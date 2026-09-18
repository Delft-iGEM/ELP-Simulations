#!/usr/bin/env bash
# A tmux session that outlives your SSH connection, with Claude Code in it.
#
#   ./work.sh            create the session (first time) or re-attach to it
#   ./work.sh status     is it running, and on which login node?
#   ./work.sh kill       end it
#
# Detach (leave it running) with  Ctrl-b  then  d.  Closing the terminal or
# losing the connection does the same thing: the session keeps running and
# Claude keeps working. Next time, log in to the SAME login node (the session
# lives on the node it was started on — this script prints which) and run
# ./work.sh again.
#
# tmux itself is not installed on DelftBlue; it lives in the conda env "tmux"
# (~/miniforge3/envs/tmux, from conda-forge), linked into ~/.local/bin.

set -euo pipefail

SESSION="claude"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NODE_FILE="$HOME/.claude-work-node"

export PATH="$HOME/.local/bin:$HOME/miniforge3/envs/tmux/bin:$PATH"
if ! command -v tmux >/dev/null 2>&1; then
    echo "tmux is missing. Install it once with:"
    echo "  ~/miniforge3/bin/conda create -y -n tmux -c conda-forge tmux"
    echo "  ln -sf ~/miniforge3/envs/tmux/bin/tmux ~/.local/bin/tmux"
    exit 1
fi

case "${1:-attach}" in
    status)
        if tmux has-session -t "$SESSION" 2>/dev/null; then
            echo "session '$SESSION' is running on $(hostname)"
            tmux list-windows -t "$SESSION"
        else
            echo "no session '$SESSION' on $(hostname)"
            [ -f "$NODE_FILE" ] && echo "last started on: $(cat "$NODE_FILE") — log in there and retry"
        fi
        ;;
    kill)
        tmux kill-session -t "$SESSION" && echo "session '$SESSION' ended"
        ;;
    attach|"")
        if ! tmux has-session -t "$SESSION" 2>/dev/null; then
            if [ -f "$NODE_FILE" ] && [ "$(cat "$NODE_FILE")" != "$(hostname)" ]; then
                echo "!  a session may still be running on $(cat "$NODE_FILE"), not on $(hostname)."
                echo "   ssh there to re-attach, or press Enter to start a new one here."
                read -r _
            fi
            hostname > "$NODE_FILE"
            tmux new-session -d -s "$SESSION" -c "$REPO" -n claude
            # A shell with Claude started in it: when Claude exits you land back in
            # the shell instead of the window closing.
            tmux send-keys -t "$SESSION:claude" "cd '$REPO' && claude" Enter
            echo "started session '$SESSION' on $(hostname) — detach with Ctrl-b d"
        else
            echo "re-attaching to '$SESSION' on $(hostname) — detach with Ctrl-b d"
        fi
        exec tmux attach-session -t "$SESSION"
        ;;
    *)
        echo "usage: $0 [attach|status|kill]"; exit 2 ;;
esac
