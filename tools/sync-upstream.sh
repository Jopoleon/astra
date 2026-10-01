#!/usr/bin/env bash
# Подтянуть обновления ASTRA из оригинала MPCDF в astra-src/ (git subtree).
#
#   tools/sync-upstream.sh
#
# Голый `git pull` с MPCDF не использовать: он притащит неочищенную историю (~1 ГБ,
# файлы по 215 МБ), у очищенной истории другие хэши. Скрипт делает так:
#   1. fetch MPCDF в полный клон .upstream/mpcdf (в .gitignore);
#   2. стоп, если MPCDF переписал свою историю (прошлая синхронизация — не предок новой);
#   3. список того, что фильтр вырежет из новых коммитов (regressions/, файлы > 5 МБ);
#   4. очистка истории тем же фильтром -> .upstream/filtered (результат детерминирован);
#   5. стоп, если очищенная история не продолжает ветку astra-upstream (хэши «поехали»);
#   6. ветка astra-upstream продвигается вперёд, фиксируется tools/upstream-sync.txt;
#   7. git subtree merge в astra-src/ текущей ветки — наши правки сохраняются.
# Ничего не пушит.
set -euo pipefail

URL="https://gitlab.mpcdf.mpg.de/git/astra.git"
BRANCH="main"
FILTER_ARGS=(--path regressions/ --invert-paths --strip-blobs-bigger-than 5M)
BIG=5000000

ROOT="$(git rev-parse --show-toplevel)"
UP="$ROOT/.upstream"
STATE="$ROOT/tools/upstream-sync.txt"
cd "$ROOT"

die() { echo "error: $*" >&2; exit 1; }
say() { echo "==> $*"; }

command -v git-filter-repo >/dev/null || die "нужен git-filter-repo: sudo apt install git-filter-repo"
[ -z "$(git status --porcelain)" ] || die "рабочее дерево не чистое — закоммитьте или уберите изменения"
[ -f "$STATE" ] || die "нет $STATE"
OLD_MPCDF="$(sed -n 's/^mpcdf=//p' "$STATE")"
OLD_FILTERED="$(sed -n 's/^filtered=//p' "$STATE")"

# 1. полный клон MPCDF
if [ ! -d "$UP/mpcdf/.git" ]; then
  say "клонирую $URL в .upstream/mpcdf (~1 ГБ, один раз)"
  mkdir -p "$UP"
  git -c core.autocrlf=false clone -q -o upstream -b "$BRANCH" "$URL" "$UP/mpcdf"
  git -C "$UP/mpcdf" config core.autocrlf false
fi
say "fetch $URL"
git -C "$UP/mpcdf" fetch -q upstream
NEW_MPCDF="$(git -C "$UP/mpcdf" rev-parse "upstream/$BRANCH")"

if [ "$NEW_MPCDF" = "$OLD_MPCDF" ]; then
  say "нового нет: MPCDF $BRANCH = ${NEW_MPCDF:0:8}"
  exit 0
fi

# 2. MPCDF не переписывал историю?
git -C "$UP/mpcdf" cat-file -e "$OLD_MPCDF^{commit}" 2>/dev/null \
  || die "прошлого синхронизированного коммита MPCDF ${OLD_MPCDF:0:8} нет в клоне — MPCDF мог переписать историю. Ничего не изменено."
git -C "$UP/mpcdf" merge-base --is-ancestor "$OLD_MPCDF" "$NEW_MPCDF" \
  || die "MPCDF переписал историю: ${OLD_MPCDF:0:8} не предок ${NEW_MPCDF:0:8}. Ничего не изменено, нужен ручной разбор."

say "новые коммиты MPCDF ${OLD_MPCDF:0:8}..${NEW_MPCDF:0:8}:"
git -C "$UP/mpcdf" log --oneline "$OLD_MPCDF..$NEW_MPCDF" | sed 's/^/    /'

# 3. что фильтр вырежет из новых коммитов
STRIPPED="$(git -C "$UP/mpcdf" rev-list --objects "$OLD_MPCDF..$NEW_MPCDF" \
  | git -C "$UP/mpcdf" cat-file --batch-check='%(objecttype) %(objectsize) %(rest)' \
  | awk -v big="$BIG" '$1=="blob" && ($3 ~ /^regressions\// || $2>big) {printf "    %.1f MB  %s\n", $2/1e6, $3}' \
  | sort -u -k3)"
if [ -n "$STRIPPED" ]; then
  echo "!!! фильтр НЕ возьмёт в наш репозиторий эти файлы из новых коммитов:"
  echo "$STRIPPED"
  echo "!!! если что-то из этого нужно — взять вручную из .upstream/mpcdf"
fi

# 4. очищенная история
git -C "$UP/mpcdf" switch -q "$BRANCH"
git -C "$UP/mpcdf" merge -q --ff-only "upstream/$BRANCH" || say "предупреждение: рабочая копия .upstream/mpcdf не обновлена (локальные изменения?)"
say "очищаю историю (filter-repo)"
rm -rf "$UP/filtered"
git -c core.autocrlf=false clone -q --no-local --single-branch -b "$BRANCH" "$UP/mpcdf" "$UP/filtered"
git -C "$UP/filtered" filter-repo --quiet --force "${FILTER_ARGS[@]}" >/dev/null
NEW_FILTERED="$(git -C "$UP/filtered" rev-parse HEAD)"

# 5. очищенная история продолжает нашу?
git fetch -q "$UP/filtered" "$BRANCH"
git merge-base --is-ancestor "$OLD_FILTERED" "$NEW_FILTERED" \
  || die "очищенная история не продолжает astra-upstream (${OLD_FILTERED:0:8} не предок ${NEW_FILTERED:0:8}): фильтр дал другие хэши. Ничего не изменено."

# 6. продвинуть astra-upstream и записать состояние
git branch -f astra-upstream "$NEW_FILTERED"
printf 'mpcdf=%s\nfiltered=%s\n' "$NEW_MPCDF" "$NEW_FILTERED" > "$STATE"
git add "$STATE"
git commit -q -m "upstream: MPCDF ${OLD_MPCDF:0:8}..${NEW_MPCDF:0:8} (filtered ${NEW_FILTERED:0:8})"

# 7. влить в astra-src/
say "subtree merge в astra-src/"
if ! git subtree merge -q --prefix=astra-src astra-upstream \
     -m "Sync ASTRA from MPCDF ${OLD_MPCDF:0:8}..${NEW_MPCDF:0:8}"; then
  echo "!!! конфликты с нашими правками: разрешите их, затем git commit" >&2
  exit 1
fi

say "готово. Проверьте и запушьте: git push origin HEAD astra-upstream"
