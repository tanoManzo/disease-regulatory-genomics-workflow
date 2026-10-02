#!/usr/bin/env python3
"""Lightweight live monitor for Claude Code background agents.

Reads the JSONL task transcripts that Claude Code streams to
/tmp/claude-*/<project-slug>/<session>/tasks/*.output and renders a
terminal status board: one row per agent with step counts, the tool it
is running right now, and its latest activity.

Usage:
  agent_monitor.py               one snapshot of the newest session
  agent_monitor.py --watch [N]   refresh every N seconds (default 5)
  agent_monitor.py --dir PATH    monitor a specific tasks directory
  agent_monitor.py --log         after the table, print each agent's
                                 last few tool calls

No dependencies beyond the standard library.
"""

import argparse
import glob
import json
import os
import sys
import time
from datetime import datetime, timezone

GREEN, YELLOW, DIM, BOLD, RESET = "\033[32m", "\033[33m", "\033[2m", "\033[1m", "\033[0m"


def find_tasks_dirs():
    dirs = glob.glob("/tmp/claude-*/*/*/tasks")
    return sorted(dirs, key=os.path.getmtime, reverse=True)


def parse_transcript(path):
    """Single pass over one JSONL transcript; returns a status dict."""
    st = {
        "file": path,
        "agent": os.path.basename(path).split(".")[0][:9],
        "task": "",
        "steps": 0,
        "tools": 0,
        "last_tool": "",
        "last_tool_ts": None,
        "open_tool": None,     # tool_use seen with no tool_result yet
        "last_text": "",
        "first_ts": None,
        "last_ts": None,
        "done": False,
    }
    pending = {}  # tool_use id -> (name, brief)
    try:
        with open(path, errors="replace") as fh:
            for line in fh:
                try:
                    o = json.loads(line)
                except json.JSONDecodeError:
                    continue
                ts = o.get("timestamp")
                if ts:
                    st["last_ts"] = ts
                    st["first_ts"] = st["first_ts"] or ts
                msg = o.get("message") or {}
                content = msg.get("content")
                if isinstance(content, str):
                    if not st["task"]:
                        st["task"] = content.strip().split("\n")[0]
                    continue
                if not isinstance(content, list):
                    continue
                for b in content:
                    bt = b.get("type")
                    if bt == "tool_use":
                        st["tools"] += 1
                        brief = tool_brief(b)
                        st["last_tool"], st["last_tool_ts"] = brief, ts
                        pending[b.get("id")] = brief
                    elif bt == "tool_result":
                        pending.pop(b.get("tool_use_id"), None)
                    elif bt == "text":
                        txt = (b.get("text") or "").strip()
                        if txt:
                            st["last_text"] = txt.split("\n")[0]
                if o.get("type") == "assistant":
                    st["steps"] += 1
    except OSError as e:
        st["last_text"] = f"unreadable: {e}"
        return st
    st["open_tool"] = next(reversed(pending.values()), None) if pending else None
    # Heuristic: a final text block with no tool still pending, and the file
    # untouched for a while, means the agent has finished its turn.
    idle = time.time() - os.path.getmtime(path)
    st["done"] = not pending and st["last_text"] and idle > 90
    return st


def tool_brief(block):
    name = block.get("name", "?")
    inp = block.get("input") or {}
    for key in ("description", "file_path", "url", "prompt", "command"):
        if key in inp and isinstance(inp[key], str):
            return f"{name}: {inp[key][:60]}"
    return name


def fmt_age(ts):
    if not ts:
        return "-"
    try:
        t = datetime.fromisoformat(ts.replace("Z", "+00:00"))
        s = int((datetime.now(timezone.utc) - t).total_seconds())
        return f"{s}s" if s < 120 else f"{s // 60}m{s % 60:02d}s"
    except ValueError:
        return "-"


def render(tasks_dir, show_log):
    files = sorted(glob.glob(os.path.join(tasks_dir, "*.output")))
    rows = [parse_transcript(f) for f in files]
    now = datetime.now().strftime("%H:%M:%S")
    print(f"{BOLD}Agent monitor{RESET}  {DIM}{tasks_dir}{RESET}  {now}\n")
    if not rows:
        print("  no agent transcripts found")
        return
    hdr = f"{'agent':<10} {'status':<8} {'steps':>5} {'tools':>5} {'last event':>10}  activity"
    print(BOLD + hdr + RESET)
    for r in rows:
        if r["done"]:
            status, color = "done", DIM
            activity = r["last_text"][:70]
        elif r["open_tool"]:
            status, color = "running", GREEN
            activity = f"in tool -> {r['open_tool'][:62]}"
        else:
            status, color = "thinking", YELLOW
            activity = (r["last_tool"] or r["last_text"])[:70]
        print(
            f"{color}{r['agent']:<10} {status:<8}{RESET} "
            f"{r['steps']:>5} {r['tools']:>5} {fmt_age(r['last_ts']):>10}  {activity}"
        )
    if show_log:
        for r in rows:
            print(f"\n{BOLD}{r['agent']}{RESET}  {DIM}{r['task'][:90]}{RESET}")
            print(f"  last tool: {r['last_tool'] or '-'}")
            print(f"  last text: {r['last_text'][:110] or '-'}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dir", help="tasks directory to monitor")
    ap.add_argument("--watch", nargs="?", const=5, type=int, metavar="SEC")
    ap.add_argument("--log", action="store_true")
    args = ap.parse_args()

    tasks_dir = args.dir or (find_tasks_dirs() or [None])[0]
    if not tasks_dir or not os.path.isdir(tasks_dir):
        sys.exit("No Claude task transcript directory found under /tmp/claude-*.")

    if args.watch:
        try:
            while True:
                os.system("clear")
                render(tasks_dir, args.log)
                time.sleep(args.watch)
        except KeyboardInterrupt:
            pass
    else:
        render(tasks_dir, args.log)


if __name__ == "__main__":
    main()
