#!/bin/sh
# Shared read-only validation for proposed commits and complete pushed histories.
guard_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
guard_author_file="$guard_dir/allowed-authors"

guard_fail() {
    printf '%s\n' "Contributor guard: $*" >&2
    return 1
}

guard_author_allowed() {
    awk -v email="$1" '
        { sub(/\r$/, ""); if (tolower($0) == tolower(email)) found = 1 }
        END { exit !found }
    ' "$guard_author_file"
}

guard_committer_allowed() {
    case "$1" in noreply@github.com) return 0 ;; esac
    guard_author_allowed "$1"
}

guard_message() {
    awk 'tolower($0) ~ /^[[:space:]]*co-authored-by[[:space:]]*:/ { print }' |
    while IFS= read -r guard_trailer; do
        guard_lower=$(printf '%s' "$guard_trailer" | tr '[:upper:]' '[:lower:]')
        case "$guard_lower" in
            *claude*|*anthropic*|*openai*|*chatgpt*|*copilot*|*codex*|*gemini*)
                guard_fail 'AI coauthor attribution is not permitted.'
                return 1 ;;
        esac
        guard_coauthor=$(printf '%s' "$guard_trailer" | sed -n 's/.*<\([^<>]*\)>[[:space:]]*$/\1/p')
        if [ -z "$guard_coauthor" ] || ! guard_author_allowed "$guard_coauthor"; then
            guard_fail "Unapproved coauthor trailer: $guard_trailer"
            return 1
        fi
    done
}

guard_commit() {
    guard_commit_oid=$1
    guard_metadata=$(git show -s --format='%ae%n%ce%n%B' "$guard_commit_oid") || return 1
    guard_author=$(printf '%s\n' "$guard_metadata" | sed -n '1p')
    guard_committer=$(printf '%s\n' "$guard_metadata" | sed -n '2p')
    if ! guard_author_allowed "$guard_author"; then
        guard_fail "Commit $guard_commit_oid has unapproved author $guard_author."
        return 1
    fi
    if ! guard_committer_allowed "$guard_committer"; then
        guard_fail "Commit $guard_commit_oid has unapproved committer $guard_committer."
        return 1
    fi
    printf '%s\n' "$guard_metadata" | sed '1,2d' | guard_message || return 1
}

guard_history() {
    guard_history_oids=$(git rev-list "$@") || return 1
    for guard_history_oid in $guard_history_oids; do
        guard_commit "$guard_history_oid" || return 1
    done
}

guard_pending_identity() {
    guard_pending_author=$(git var GIT_AUTHOR_IDENT | sed -n 's/.*<\([^<>]*\)>.*/\1/p') || return 1
    guard_pending_committer=$(git var GIT_COMMITTER_IDENT | sed -n 's/.*<\([^<>]*\)>.*/\1/p') || return 1
    if ! guard_author_allowed "$guard_pending_author"; then
        guard_fail "Use an approved author email, not $guard_pending_author."
        return 1
    fi
    if ! guard_committer_allowed "$guard_pending_committer"; then
        guard_fail "Use an approved committer email, not $guard_pending_committer."
        return 1
    fi
}

guard_pending_history() {
    if git rev-parse --verify HEAD >/dev/null 2>&1; then
        guard_history HEAD || return 1
    fi
    guard_merge_file=$(git rev-parse --git-path MERGE_HEAD) || return 1
    if [ -f "$guard_merge_file" ]; then
        while IFS= read -r guard_merge_oid; do
            guard_history "$guard_merge_oid" || return 1
        done < "$guard_merge_file"
    fi
}
