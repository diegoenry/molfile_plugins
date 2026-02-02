#!/bin/bash
# Launcher script for GROMACS H5MD Browser Web Server
# Usage: ./run_server.sh [directory] [--port PORT]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARENT_DIR="$(dirname "$SCRIPT_DIR")"
VENV_DIR="$PARENT_DIR/venv"

# Check if venv exists, if not try to create it
if [ ! -d "$VENV_DIR" ]; then
    echo "Virtual environment not found. Creating one..."
    python3 -m venv "$VENV_DIR"
    "$VENV_DIR/bin/pip" install --upgrade pip
    "$VENV_DIR/bin/pip" install -r "$SCRIPT_DIR/requirements.txt"
fi

# Install flask if not present
"$VENV_DIR/bin/pip" show flask > /dev/null 2>&1 || "$VENV_DIR/bin/pip" install flask

# Use parent directory as default (where h5md files are)
DIR="${1:-$PARENT_DIR}"

echo "Starting H5MD Browser Web Server..."
echo "  Directory: $DIR"
echo "  URL: http://127.0.0.1:5000"
echo ""
echo "Press Ctrl+C to stop the server"
echo ""

exec "$VENV_DIR/bin/python" "$SCRIPT_DIR/app.py" "$DIR" "${@:2}"
