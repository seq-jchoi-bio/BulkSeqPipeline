#!/bin/bash

expand_input_path() {
    local input="${1:-}"
    local user rest home_dir

    case "$input" in
        \~)
            printf "%s\n" "$HOME"
            ;;
        \~/*)
            printf "%s/%s\n" "$HOME" "${input:2}"
            ;;
        \~*/*)
            user="${input%%/*}"
            user="${user#~}"
            rest="${input#*/}"
            if [[ "$user" =~ ^[A-Za-z0-9._-]+$ ]]; then
                if command -v getent >/dev/null 2>&1; then
                    home_dir="$(getent passwd "$user" | cut -d: -f6)"
                elif command -v dscl >/dev/null 2>&1; then
                    home_dir="$(dscl . -read "/Users/$user" NFSHomeDirectory 2>/dev/null | awk '{print $2}')"
                else
                    home_dir=""
                fi

                if [ -n "$home_dir" ]; then
                    printf "%s/%s\n" "$home_dir" "$rest"
                else
                    printf "%s\n" "$input"
                fi
            else
                printf "%s\n" "$input"
            fi
            ;;
        *)
            printf "%s\n" "$input"
            ;;
    esac
}

detect_cpu_total() {
    if command -v nproc >/dev/null 2>&1; then
        nproc
    else
        getconf _NPROCESSORS_ONLN
    fi
}

detect_memory_mb() {
    local mem_bytes
    local mem_value mem_unit

    if command -v free >/dev/null 2>&1; then
        free -m | awk '/^Mem:/ {print $2; exit}'
    elif command -v sysctl >/dev/null 2>&1; then
        mem_bytes="$(sysctl -n hw.memsize 2>/dev/null || echo 0)"
        if [[ "$mem_bytes" =~ ^[0-9]+$ ]] && [ "$mem_bytes" -gt 0 ]; then
            echo $(( mem_bytes / 1024 / 1024 ))
        elif [ -r /proc/meminfo ]; then
            awk '/^MemTotal:/ {print int($2 / 1024); exit}' /proc/meminfo
        elif command -v system_profiler >/dev/null 2>&1; then
            read -r mem_value mem_unit < <(system_profiler SPHardwareDataType 2>/dev/null | awk -F': ' '/Memory:/ {print $2; exit}')
            if [[ "$mem_value" =~ ^[0-9]+$ ]]; then
                case "$mem_unit" in
                    TB) echo $(( mem_value * 1024 * 1024 )) ;;
                    GB) echo $(( mem_value * 1024 )) ;;
                    MB) echo "$mem_value" ;;
                    *) echo 0 ;;
                esac
            else
                echo 0
            fi
        else
            echo 0
        fi
    elif [ -r /proc/meminfo ]; then
        awk '/^MemTotal:/ {print int($2 / 1024); exit}' /proc/meminfo
    else
        echo 0
    fi
}

format_memory_gb() {
    local mem_mb="${1:-0}"
    awk -v mb="$mem_mb" 'BEGIN { printf "%.1f", mb / 1024 }'
}

is_positive_integer() {
    [[ "${1:-}" =~ ^[0-9]+$ ]] && [ "$1" -ge 1 ]
}

is_valid_samtools_mem() {
    [[ "${1:-}" =~ ^[0-9]+[KkMmGg]?$ ]]
}

optimize_resources() {
    local requested_cpu="${1:?missing requested CPU}"
    local cpu_total="${2:?missing total CPU}"
    local mem_mb="${3:-0}"
    local reserve_mb usable_mb max_cpu_by_mem cpu sort_threads sort_mem_mb sort_budget_mb

    cpu="$requested_cpu"
    if [ "$cpu" -gt "$cpu_total" ]; then
        cpu="$cpu_total"
    fi

    sort_threads="$cpu"
    sort_mem_mb=1024

    if [ "$mem_mb" -gt 0 ]; then
        reserve_mb=$(( mem_mb / 4 ))
        if [ "$reserve_mb" -lt 2048 ]; then
            reserve_mb=2048
        fi
        if [ "$reserve_mb" -ge "$mem_mb" ]; then
            reserve_mb=$(( mem_mb / 3 ))
        fi

        usable_mb=$(( mem_mb - reserve_mb ))
        if [ "$usable_mb" -lt 1024 ]; then
            usable_mb=1024
        fi

        max_cpu_by_mem=$(( usable_mb / 2048 ))
        if [ "$max_cpu_by_mem" -lt 1 ]; then
            max_cpu_by_mem=1
        fi
        if [ "$cpu" -gt "$max_cpu_by_mem" ]; then
            cpu="$max_cpu_by_mem"
        fi

        sort_threads="$cpu"
        sort_budget_mb=$(( usable_mb * 70 / 100 ))
        sort_mem_mb=$(( sort_budget_mb / sort_threads ))
        if [ "$sort_mem_mb" -lt 512 ]; then
            sort_mem_mb=512
        elif [ "$sort_mem_mb" -gt 2048 ]; then
            sort_mem_mb=2048
        fi
    fi

    OPT_CPU="$cpu"
    OPT_SORT_THREADS="$sort_threads"
    OPT_SORT_MEM="${sort_mem_mb}M"
}

prepare_resource_settings() {
    local requested_cpu="${1:-2}"
    local context="${2:-Pipeline step}"
    local provided_sort_mem="${3:-}"
    local cpu_total mem_total_mb

    if ! is_positive_integer "$requested_cpu"; then
        echo "Warning: Invalid thread value '$requested_cpu'. Using 2 threads before system safety checks."
        requested_cpu=2
    fi

    cpu_total="$(detect_cpu_total)"
    mem_total_mb="$(detect_memory_mb)"
    optimize_resources "$requested_cpu" "$cpu_total" "$mem_total_mb"

    RESOURCE_REQUESTED_CPU="$requested_cpu"
    RESOURCE_CPU_TOTAL="$cpu_total"
    RESOURCE_MEM_TOTAL_MB="$mem_total_mb"
    RESOURCE_CPU="$OPT_CPU"
    RESOURCE_SORT_THREADS="$OPT_SORT_THREADS"
    if [ -n "$provided_sort_mem" ]; then
        if is_valid_samtools_mem "$provided_sort_mem"; then
            RESOURCE_SORT_MEM="$provided_sort_mem"
        else
            echo "Warning: Invalid samtools sort memory value '$provided_sort_mem'. Using calculated value '$OPT_SORT_MEM'."
            RESOURCE_SORT_MEM="$OPT_SORT_MEM"
        fi
    else
        RESOURCE_SORT_MEM="$OPT_SORT_MEM"
    fi

    if [ "$RESOURCE_REQUESTED_CPU" -ne "$RESOURCE_CPU" ]; then
        echo "==========================================="
        echo "   Resource adjustment: ${context}"
        echo "==========================================="
        echo "   Requested threads: ${RESOURCE_REQUESTED_CPU}"
        echo "   Available threads: ${RESOURCE_CPU_TOTAL}"
        if [ "$RESOURCE_MEM_TOTAL_MB" -gt 0 ]; then
            echo "   Total memory     : ~$(format_memory_gb "$RESOURCE_MEM_TOTAL_MB") GB"
        fi
        echo "   Using threads    : ${RESOURCE_CPU}"
        echo "==========================================="
    fi
}
