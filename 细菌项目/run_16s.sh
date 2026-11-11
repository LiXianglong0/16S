#!/bin/bash
PROJECT_DIR=$(cd "$(dirname "$0")" && pwd)
bash "$PROJECT_DIR/../pipeline/run_16s.sh" "$PROJECT_DIR" "$@"
