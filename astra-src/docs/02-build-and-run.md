# 02. Сборка и запуск ASTRA 8: конвейер от `as_exe` до CDF

**О чём глава.** Как именно устроены установка, сборка и запуск ASTRA 8: что делает `install.sh`; как выбирается площадка; каждый файл `repo_user_area/exe/` (`as_exe`, `as_exe.py` со всеми опциями, `Build`, `Makefile`, `Makeastra`, `Makesub`, `Makexpr`, `astra_rc*`, `greenMatrices.py`, `version`, `wr_nml`, `README`); внешние модули и раскладка `$ASTRA_EXT` (`modules_install.sh`, `modules_inst_loc.sh`); `requirements.txt`; какие каталоги появляются при сборке и расчёте; полная последовательность одного запуска; batch, GUI и SLURM; переменные окружения; пересборка и очистка; типичные ошибки. Всё сверено с кодом `astra-src` (`67281112`, 2026-10-01) и с рабочей установкой `~/astra/a8` (`e66f6ba1` + патчи установщика).

Что **не** повторяется здесь:

- проблемы сборки на Ubuntu и как их решает наш установщик `install_astra8.sh` — [`../../installer/docs/DETAILED.md`](../../installer/docs/DETAILED.md) (разделы 3, 4, 8);
- разбор всех `platform/env.*` и `get_platform` по площадкам — [`../../docs/astra/platforms.md`](../../docs/astra/platforms.md);
- общий обзор и карта дерева — [01-overview.md](01-overview.md);
- что генерирует парсер и как устроен язык модели — [03-model-language.md](03-model-language.md);
- внешние физические модули по существу — [07-heating-and-external-modules.md](07-heating-and-external-modules.md);
- вывод, GUI, регрессии — [08-graphics-and-output.md](08-graphics-and-output.md).

## Содержание

1. [Конвейер одним взглядом](#1-конвейер-одним-взглядом)
2. [`install.sh`](#2-installsh)
3. [Выбор площадки: `get_platform` и `platform/env.*`](#3-выбор-площадки-get_platform-и-platformenv)
4. [`repo_user_area/exe/`: файл за файлом](#4-repo_user_areaexe-файл-за-файлом)
5. [Внешние модули и `$ASTRA_EXT`](#5-внешние-модули-и-astra_ext)
6. [`requirements.txt` и реальные зависимости](#6-requirementstxt-и-реальные-зависимости)
7. [Каталоги и файлы, которые появляются при работе](#7-каталоги-и-файлы-которые-появляются-при-работе)
8. [Полная последовательность одного запуска](#8-полная-последовательность-одного-запуска)
9. [Batch, GUI, SLURM, flight simulator](#9-batch-gui-slurm-flight-simulator)
10. [Пересборка и очистка](#10-пересборка-и-очистка)
11. [Переменные окружения](#11-переменные-окружения)
12. [Типичные ошибки](#12-типичные-ошибки)
13. [Где смотреть в исходниках](#13-где-смотреть-в-исходниках)

---

## 1. Конвейер одним взглядом

```
exe/as_exe [опции]                               bash: площадка, env, PATH к python
 └─ python3 exe/as_exe.py [опции]
     ├─ greenMatrices.main()                     exp/cnf/*_description_in.json → *_full.json
     ├─ пишет tmp/astra.nml (&astra_log)         параметры запуска; источник «последнего расчёта»
     ├─ os.system("exe/Build <m> <e> <RUN|BGD> <fs> <nodes> <workflow> [-W]")
     │   ├─ cd $AWD; get_platform; source platform/env.<pl>; source exe/astra_rc[_<comp>]
     │   ├─ ulimit -s 100000
     │   ├─ make -f exe/Makefile MODEL=<m> FEXP=<e> ATASK=bin/<m>_{gui|batch}.exe WORKFLOW=...
     │   │   ├─ make -f exe/Makeastra ... sub={misc,astra,feqis,graph,nbi} modules   (.mod)
     │   │   ├─ make -f exe/Makeastra ... sub=... lib     → lib/lib{misc,astra,feqis,graph,nbi}.a
     │   │   ├─ make -f exe/Makesub sub=sbr | sub=fnc     → lib/libsbr.a, lib/libfnc.a
     │   │   ├─ pyparse/parser_main.py -equ equ/<m> -exp exp/<e>   → src/tmp/*.f90
     │   │   └─ $(FC) astra_main.F90 + stepup.F90 + src/tmp/*.f90 + libs → bin/<m>_*.exe
     │   ├─ make -f exe/Makexpr                  xpr/tglfi, xpr/qlki, xpr/neo
     │   ├─ проверка конфликтов имён libsbr.a / libfnc.a
     │   └─ запуск bin/<m>_*.exe  (или sbatch, если есть SLURM и -batch)
     ├─ уборка IPC: kill + ipcrm по tmp/<e><m>-{0,1,2}.ipc
     └─ postProcess/json2cdf.json_concat()       ncdf_out/<e><m>-*.json → ncdf_out/<e><m>.CDF
```

Ключевые свойства:

- исполняемый файл **свой для каждой модели** и пересобирается, если изменилась модель, `fml/`, `pyparse/`, `exe/astra_rc`, `sbr/`, `fnc/` или исходники ядра (раздел 10);
- `as_exe` **всегда возвращает 0** (результат `os.system` в `as_exe.py:131` игнорируется) — успех определяется по строке `>>> ASTRA normal exit >>>` и наличию CDF ([`DETAILED.md`](../../installer/docs/DETAILED.md), разд. 3.9);
- всё выполняется из корня установки `$AWD`: исполняемый файл открывает `tmp/astra.nml`, `equ/log/...`, `exp/...` по относительным путям (`src/astra/read_input.f90`, строки 94–152).

## 2. `install.sh`

Файл в корне, 88 строк. Это **upstream-установщик для площадок IPP и партнёров**, не для чистой Ubuntu (там он падает — [`DETAILED.md`](../../installer/docs/DETAILED.md), разд. 3). Аргумент один — необязательный `-safe`.

Пошагово (`install.sh`):

| строки | действие |
|---|---|
| 3–4 | `AWD` = абсолютный путь каталога скрипта, экспортируется |
| 6–7 | `chmod 744 get_platform`; `platform=$(get_platform)` |
| 13–29 | цикл по **всем** подкаталогам `repo_user_area/` (`equ`, `exe`, `exp`, `fml`, `fnc`, `preProcess`, `pyparse`, `sbr`, `udb`, `xpr`): `cp -r repo_user_area/$DIR .` |
| 15 | условие копирования без вопроса: каталога ещё нет, **или** режим не `-safe`, **или** это `equ`, `exp`, `udb` |
| 18–27 | иначе (только `-safe`, каталог уже есть): вопрос `'$DIR': backup + copy from repo_user_area_dir? y/n`; при `y` — `rm -rf ${DIR}_backup`, `cp -r $DIR ${DIR}_backup`, затем копирование |
| 31–63 | только `-safe`: вопросы про RABBIT, TORBEAM, QuaLiKiZ, QLKNN, NEO, TGLF. Ответ `n`: в `equ/flux_feqis` `RABBIT`→`!RABBIT`, `TORBA`→`!TORBA` (`sed`); в `exe/Makexpr` из строки `all: directories` вырезаются `$(XPR)/qlki`, `$(XPR)/neo`, `$(XPR)/tglfi`. Ответ про QLKNN читается, но ни на что не влияет |
| 66–67 | `source platform/env.$platform`; печать `$ASTRA_EXT` |
| 69 | `mkdir -p tmp` |
| 71–76 | `chmod 744` для `exe/Build`, `exe/as_exe`, `exe/wr_nml`, `pyparse/parser_main.py`, `clean.sh`, `compReg.sh` (в git `exe/*` хранятся с режимом `100644`, поэтому без этого шага `exe/as_exe` не запускается) |
| 78–87 | `make -f exe/Makefile clean` (в `-safe` — после вопроса) |
| 88 | тестовый расчёт `exe/as_exe -m flux_feqis -v aug34954 -s 4 -e 5` (GUI, т.к. без `-batch`) |

Что важно знать:

- **Перезапись файлов пользователя.** `cp -r` перезаписывает одноимённые файлы. Без `-safe` перезаписываются все каталоги, включая изменённые пользователем `sbr/`, `fml/`, `exe/`. Даже в `-safe` каталоги `equ`, `exp`, `udb` копируются **без вопроса** (строка 15), т.е. изменённые поставляемые примеры (например, своя правка `equ/flux_feqis`) будут затёрты. Свои файлы с другими именами не удаляются.
- **`pyparse` тоже копируется**: `repo_user_area/pyparse/iondens_ass.py` попадает в репозиторный `pyparse/` (файл игнорируется git).
- Опечатка в сообщении (строка 26): `repo_use_area`.
- Тестовый расчёт требует RABBIT и TORBEAM, если не ответить `n` в `-safe` ([`DETAILED.md`](../../installer/docs/DETAILED.md), разд. 3.7), и окна X11.
- Внешние модули `install.sh` не ставит — для этого `modules_install.sh` / `modules_inst_loc.sh` (раздел 5).
- Наш установщик `install_astra8.sh` не вызывает `install.sh`, а повторяет его шаги идемпотентно и с поправками ([`DETAILED.md`](../../installer/docs/DETAILED.md), разд. 4).

## 3. Выбор площадки: `get_platform` и `platform/env.*`

`get_platform` (77 строк) печатает имя площадки: сначала по пути установки (`/localhome*` → `goblin`, `/home/ITER*` → `iter`, `/home/matlab*` → `docker`, `/global/homes*` → `perlmutter`, `/compass*` → `cz`), затем по `uname -n` (`*mit.edu`/`node*` → `mit`, `omega*`, `tok*`, `viz*` → `gway`, `iris*`, ...), затем `Darwin` → `darwin`, иначе `mypc`. Если задана переменная `ASTRA_COMPILER`, к имени добавляется `_${ASTRA_COMPILER}` (строки 72–75). Подробный разбор правил, мёртвых записей (`goblin`, `iris` — файлов `env.goblin`, `env.iris` нет) и каждого `env.*` — в [`platforms.md`](../../docs/astra/platforms.md).

Кто вызывает `get_platform` и читает `platform/env.<pl>`:

| кто | зачем |
|---|---|
| `install.sh` | `ASTRA_EXT` для печати; переменные для `make clean` |
| `exe/as_exe` | `PYTHON_BINDIR` → `PATH`; остальное наследует `as_exe.py` и `Build` |
| `exe/Build` | компилятор, флаги, `MKL_LIBDIR`, `ASTRA_EXT`, `SL_QUEUE` |
| `exe/Makefile`, `exe/Makexpr` | только имя (`PLATFORM := $(shell $(AWD)/get_platform)`), для ветки `darwin` |
| `exe/astra_rc` | только имя, для ветки `darwin` |
| `compReg.sh`, `modules_install.sh`, `modules_inst_loc.sh` | то же, что `as_exe` / установка модулей |

Контракт переменных, которые должен задать `env.<pl>` (минимум для сборки): `FC`, `CC`, `FC_FLAGS`, `LD_FLAGS`, `ASTRA_EXT`, `MKL_LIBDIR` (кроме `darwin`, где `OPENBLAS_LIBDIR`); для IPC-программ — `MPIFC`, `MPICC`, а при Intel — `COMP_LIBDIR`; желательно `PYTHON_BINDIR`. Для Ubuntu это файл `platform/env.ubuntu`, который пишет установщик (текст — [`DETAILED.md`](../../installer/docs/DETAILED.md), разд. 3.2).

Установщик дополнительно вставляет в начало `get_platform` переключатель `ASTRA_PLATFORM` (`~/astra/a8/get_platform`, строки 3–4 после патча): при заданной переменной скрипт сразу печатает её и выходит. Побочный эффект: суффикс `_${ASTRA_COMPILER}` при этом не добавляется, т.е. опция `-c` у `as_exe` на установке через `astra8` игнорируется при выборе `env.*` (но `Build` всё равно возьмёт `exe/astra_rc_${ASTRA_COMPILER}`, см. 4.3).

## 4. `repo_user_area/exe/`: файл за файлом

После установки каталог копируется в `$AWD/exe/` и дальше используется **копия**. 14 файлов:

| файл | тип | роль |
|---|---|---|
| `as_exe` | bash | точка входа: площадка, `PATH`, вызов `as_exe.py` |
| `as_exe.py` | Python | разбор опций, `tmp/astra.nml`, вызов `Build`, уборка IPC, JSON → CDF |
| `Build` | bash | окружение, `make`, проверка конфликтов, запуск или `sbatch` |
| `Makefile` | GNU make | главная сборка: библиотеки, парсер, линковка |
| `Makeastra` | GNU make | компиляция одного каталога `src/<sub>` в модули и `lib<sub>.a` |
| `Makesub` | GNU make | компиляция `sbr/` или `fnc/` в `libsbr.a` / `libfnc.a` |
| `Makexpr` | GNU make | IPC-программы `xpr/tglfi`, `xpr/qlki`, `xpr/neo` |
| `astra_rc` | bash | пути к внешним модулям в `$ASTRA_EXT` (по умолчанию) |
| `astra_rc_gcc`, `astra_rc_ifx` | bash | то же для `ASTRA_COMPILER=gcc` / `ifx` |
| `greenMatrices.py` | Python | матрицы Грина катушек для FEQIS |
| `version` | текст | баннер с номером версии |
| `wr_nml` | bash | запуск редактора namelist из TRVIEW (только IPP) |
| `README` | текст | 10 строк о каталогах `mod/`, `lib/` |

### 4.1. `as_exe`

13 строк (`repo_user_area/exe/as_exe`):

1. `AWD` = родитель каталога скрипта (не экспортируется);
2. `platform=$(get_platform)`; `source platform/env.$platform`;
3. если `PYTHON_BINDIR` нет в `PATH` — добавить его в начало;
4. `python3 ${AWD}/exe/as_exe.py "$@"`.

Код возврата скрипта — код возврата `as_exe.py`, то есть 0 почти всегда.

### 4.2. `as_exe.py`

152 строки. Импортирует `numpy` и `greenMatrices` (лежит рядом, поэтому импорт работает при запуске из любого каталога), затем `postProcess.json2cdf`.

**Порядок действий** (`__main__`, строки 38–152):

1. `greenMatrices.main()` (строка 43) — пересчёт описаний машин в `exp/cnf/` (раздел 4.9). Выполняется **при каждом** вызове, даже для моделей без FEQIS.
2. `parse_alog()` (строки 24–35) читает `tmp/astra.nml` построчно, делит по `=`, снимает кавычки. Если файла нет — `None`.
3. Строит `argparse` (строки 45–68). Если `tmp/astra.nml` нет — `-m` и `-v` **обязательны**; если есть — значения по умолчанию берутся из него.
4. Определяет `rtype`: `BGD` при `-batch`, иначе `RUN`.
5. Пишет `tmp/astra.nml` (строки 101–121), namelist `&astra_log`.
6. Вызывает `exe/Build <equ> <exp> <rtype> <fs> <nodes> <fpath> [-W]` через `os.system` (строки 127–131).
7. Убирает «зависшие» IPC-процессы (строки 133–147): для `ipcId` 0, 1, 2 читает `tmp/<exp><equ>-<ipcId>.ipc`, берёт идентификаторы разделяемой памяти из строк 7–9 и таблицу PID начиная с 11-й строки; `kill -9` каждому PID, `ipcrm -m` каждому сегменту.
8. `j2nc.json_concat(exp+equ)` (строки 148–152) в `try/except: pass` — любая ошибка конвертации **молча** проглатывается.

**Все опции командной строки**:

| опция | длинная | тип, по умолчанию (без `tmp/astra.nml`) | по умолчанию (есть `tmp/astra.nml`) | запоминается | куда идёт |
|---|---|---|---|---|---|
| `-m` | `--equ` | строка, **обязательна** | `equ_file` | да | имя модели `equ/<m>`; `Build` $1 |
| `-v` | `--exp` | строка, **обязательна** | `exp_file` | да | имя эксперимента `exp/<e>`; `Build` $2 |
| `-s` | `--tbeg` | float, `0.1` | `tbeg_nml` | да | `tbeg_nml` → `TSTART` (перекрывает `equ/log`) |
| `-e` | `--tend` | float, `10.` | `tend_nml` | да | `tend_nml` → `TEND` |
| `-dev` | `--DeviceName` | строка, `aug` | `machine` из файла, иначе `aug` | да | `machine`: `exp/cnf/<dev>_description_in.json`, запасные `exp/nml/<dev>`, `exp/nbi/<dev>`, файл FEQIS `exp/cnf/<dev>_description_full.json` |
| `-resize` | `--resize` | float, `1.` | `resize` из файла, иначе 1 | да | масштаб окна GUI |
| `-nodes` | `--nodes` | int, `1` | `1` | нет | `Build` $5 → `N_NODES` (используется в `src/astra/a2vmec.f90:313`) |
| `-batch` | — | флаг | — | **нет** | `TASK="BGD"`, сборка `bin/<m>_batch.exe` |
| `-tpause` | `--tpause` | float, `1.e4` | `1.e4` | нет | `tpause_nml` → `TPAUSE` (пауза GUI в момент времени) |
| `-debug` | — | флаг | — | нет | `debug = 1` (модуль `src/astra/debugger.f90`) |
| `-re` | `--restart` | int, `0` | `0` | нет | `restart = N`: начальное состояние из `ncdf_out/<e><m>-N.json` (`src/astra/astra_main.F90:63`, `src/astra/json_rw.f90:737`) |
| `-fs` | — | флаг | — | **нет** | `flightsim = .True.`; `Build` $4; принудительно batch-сборка |
| `-workflow` | `--fpath` | путь, `$AWD/src/astra/stepup.F90` | то же | нет | файл, компилируемый **вместо** `stepup.F90` (`Makefile`, `ASTRA_SRC`) |
| `-W` | `--waitslurm` | флаг | — | нет | вместе с `-batch` — `sbatch -W` (ждать завершения задания) |
| `-c` | `--compiler` | строка (`gcc`, `ifx`) | `compiler` из файла | да | `ASTRA_COMPILER` в окружение → `get_platform` даёт `<pl>_<comp>`, `Build` читает `exe/astra_rc_<comp>` |
| `-h` | `--help` | — | — | — | справка argparse |

Сокращения argparse работают: `-b` = `-batch`, `-f` = `-fs` (проверено на Python 3.12 с тем же набором опций); `-d` неоднозначен (`-dev`/`-debug`).

**Пример `tmp/astra.nml`** (`~/astra/a8/tmp/astra.nml` после теста установщика):

```
&astra_log

exp_file  = "aug34954"
equ_file  = "flux_feqis"
tbeg_nml  =   4.0000
tend_nml  =   5.0000
TASK      = "BGD"
machine   = "aug"
debug     = 0
flightsim = .False.
resize    =   1.0000
workflow  = "/home/egor/astra/a8/src/astra/stepup.F90"
restart   = 0
tpause_nml= 10000.0000

/
```

Этот же файл читает исполняемый файл (`src/astra/read_input.f90::read_nml`, namelist `astra_log` со списком `equ_file, exp_file, task, machine, debug, tbeg_nml, tend_nml, tpause_nml, resize, restart, flightsim, workflow`).

Тонкости:

- `-batch` и `-fs` не запоминаются: `astra8` без аргументов после batch-расчёта откроет **окно**.
- При `-c` в файл пишется `compiler = gcc` без кавычек, а в namelist `astra_log` (Fortran) переменной `compiler` нет. Чтение идёт с `iostat=ios` без проверки (`read_input.f90`, строки 135–138), а строка стоит последней, поэтому, по-видимому, все остальные значения успевают прочитаться. Не проверено запуском.
- `parse_alog` делит строку по `=`, поэтому путь `workflow` со знаком `=` сломает разбор.

### 4.3. `Build`

122 строки (`repo_user_area/exe/Build`). Аргументы: `$1` модель, `$2` эксперимент, `$3` `RUN|BGD`, `$4` `True|False` (flight simulator), `$5` число узлов, `$6` файл workflow, `$7` необязательный `-W`.

| строки | действие |
|---|---|
| 8 | `export N_NODES=$5` |
| 16–19 | `export AWD=<корень>`, `cd $AWD` |
| 22–23 | `source platform/env.$(get_platform)` |
| 24–29 | `source exe/astra_rc_${ASTRA_COMPILER}`, если переменная задана, иначе `source exe/astra_rc` |
| 31 | `ulimit -s 100000` (иначе `unlimited`) — большой стек для автоматических массивов Fortran |
| 35–42 | batch (`FSIM == True` или `TASK == BGD`): `unset X11`, `ATASK=bin/<m>_batch.exe`; иначе `export X11=y`, `ATASK=bin/<m>_gui.exe` |
| 43–49 | `make -f exe/Makefile MODEL=... FEXP=... ATASK=... WORKFLOW=...`; при ошибке `Error from exe/Makefile-compilation, exiting`, код 1 |
| 51–58 | `make -f exe/Makexpr`; при ошибке `Error from exe/Makexpr-compilation, exiting`, код 2 |
| 62–79 | конфликт имён: каждый объект из `ar -t lib/libsbr.a` сравнивается с каждым из `lib/libfnc.a`; совпадение → `Conflicting <obj>: in libsbr.a and libfnc.a`, `No ASTRA simulation started`, код 3 |
| 83–107 | `BGD`: если `sbatch` **нет** — `Running ... without graphic interactive window` и запуск; если **есть** — пишется `exe/batch.sh` (`#SBATCH -J astra_$USER`, `${SL_QUEUE}`, экспорт `ASTRA_EXT`, `PATH`, `LD_LIBRARY_PATH`, `MAX_NWORKERS=128`) и `sbatch $WAIT_SLURM exe/batch.sh` |
| 108–121 | GUI: `MAX_NWORKERS = nproc / (потоков на ядро из lscpu)`, `Max nworkers ...`, `Running ...`, запуск, затем **безусловно** `echo ASTRA normal exit` |

Наблюдения:

- В GUI-ветке строка `ASTRA normal exit` (без `>>>`) печатается скриптом при любом исходе. Признак успеха — только `>>> ASTRA normal exit >>>` из `src/astra/astra_main.F90:131`.
- В локальной batch-ветке `MAX_NWORKERS` **не задаётся**, и IPC-интерфейсы (`sbr/tglf_ipc.f90:88`, `sbr/qlk_ipc.f90:77`, `sbr/neo_ipc.f90:81`) используют своё запасное значение (`MAX_NWORKERS not set, using fallback`).
- Коды 1/2/3 до пользователя не доходят: `as_exe.py` их игнорирует.
- Если на машине установлен SLURM (`sbatch` в `PATH`), `-batch` **отправляет задание в очередь**, а не считает локально; `SL_QUEUE` берётся из `env.<pl>` (для неизвестной площадки — пустой).

### 4.4. `Makefile`

237 строк (`repo_user_area/exe/Makefile`). Требует `AWD` в окружении (строки 1–3: `$(error setenv/export AWD first)`).

**Пути и переменные** (строки 5–35): `LIB_DIR=$(AWD)/lib`, `TSK_DIR=$(AWD)/bin`, `MOD_DIR=$(AWD)/mod`, `OBJ_DIR=$(AWD)/obj`, `TMP_DIR=$(AWD)/src/tmp`, `EQU=$(AWD)/equ/$(MODEL)`, `EXP=$(AWD)/exp/$(FEXP)`. К `FC_FLAGS` и `LD_FLAGS` добавляется `-I$(AWD)` — чтобы `sbr/*.f90` и `src/tmp/*.f90` могли `include` файлы `fml/...`.

**Генерируемые файлы** (`TMP_SRC`, строки 39–46): `associate_pointers.f90`, `inivar.f90`, `detvar.f90`, `postep.f90`, `ininam.f90`, `init_converge_step.f90`, `eqns_inc.f90`; при `X11` — ещё `astra_out.f90`, `set_graph_names.f90` (строка 67). Парсер всегда пишет все 9 файлов плюс `model.txt`, `declar.fml`, `declar.fnc` (`pyparse/parser_main.py::write_tmp`), но в batch-линковку два GUI-файла не входят.

**X11** (строки 52–72): если `X11` задан — `-DX11`, `-lgraph` и путь к libX11: `/opt/X11/lib` на `darwin`, иначе `/usr/lib/$(uname -m)-linux-gnu`, если там есть `libX11.so*` (Debian/Ubuntu), иначе `/usr/lib64` (SUSE/RHEL). Комментарий в файле говорит, что раньше для x86_64 Debian ветки не было и GUI не линковался.

**Объекты ядра** (строки 76–151): списки модулей в порядке зависимостей — `AMODS` (25 модулей `src/astra`), `AOBJS` (`ipc_control.o`), `FMODS`/`FOBJS` (FEQIS), `GMODS`/`GOBJS` (графика), `MMODS` (misc), `NMODS` (NBI). Порядок в списке и есть порядок компиляции модулей — автоматического анализа зависимостей нет.

**Линковка** (строки 159–166):

```
MKL_LIB = -Wl,--start-group $(MKL_LIBDIR)/libmkl_intel_lp64.a $(MKL_LIBDIR)/libmkl_core.a \
          $(MKL_LIBDIR)/libmkl_sequential.a -Wl,--end-group      # на darwin: -L$(OPENBLAS_LIBDIR) -lopenblas
LIB = -lsbr -lnbi $(X11_LIB) $(QLKNN_LIB) $(TGLF_LIB) $(TRAVIS_LIB) $(TORBEAM_LIB) $(RABBIT_LIB) \
      $(GRAPH_FLAG) -lastra -lfnc -lfeqis -lmisc $(SPIDER_LIB) $(JSON_LIB) -L$(LIB_DIR) \
      $(MKL_LIB) $(RT_LIB) $(TORIC_LIB)
```

`-lX11` стоит раньше `-lgraph` — отсюда проблема `--as-needed` на Ubuntu ([`DETAILED.md`](../../installer/docs/DETAILED.md), разд. 3.4). `MKL_LIB` присваивается внутри Makefile, поэтому заменить его можно только через переменную командной строки make или `MAKEFLAGS` (так делает `env.ubuntu`: `MAKEFLAGS="MKL_LIB=-lopenblas"`). `TRAVIS_LIB` нигде не задаётся (остаток). Закомментированная строка 156 — вариант «FOR STELLARATOR» с `LIBSTELL` и `-lstdc++`.

**Цели** (строки 180–237):

| цель | что делает |
|---|---|
| `all` (по умолчанию) | `directories modules $(ATASK)` |
| `directories` | `mkdir -p bin src/tmp lib ncdf_out dat` |
| `modules` | 5 вызовов `make -f exe/Makeastra ... modules` для `misc`, `astra`, `feqis`, `graph`, `nbi` — компиляция модулей в нужном порядке и перенос `.mod` |
| `$(ATASK)` | зависит от `lib{graph,astra,misc,nbi,feqis}.a`, `$(SPIDER_LIB)`, `libsbr.a`, `libfnc.a`, `exe/astra_rc`, `equ/$(MODEL)`, `astra_main.F90`, `$(WORKFLOW)`, всех `pyparse/*.py`, всех `fml/*`. Рецепт: `$(PY_TMP_GEN)` (парсер), затем `$(FC) $(LD_FLAGS) -o $(ATASK) astra_main.F90 $(WORKFLOW) $(TMP_SRC) -I... $(LIB) $(X11_FLAG)` |
| `lib/libsbr.a` | от `sbr/*` → `make -f exe/Makesub sub=sbr` |
| `lib/libfnc.a` | от `fnc/*.f90` → `make -f exe/Makesub sub=fnc` |
| `lib/libmisc.a` и др. | от `src/<sub>/*` → `make -f exe/Makeastra ... sub=<sub> lib` |
| `clean` | `rm -f *.mod`; `make -f exe/Makexpr clean`; `rm -rf mod bin obj lib` |

Главная программа (`astra_main.F90`), шаг (`stepup.F90` или `-workflow`) и сгенерированный код компилируются **прямо в команде линковки**, без объектных файлов, — поэтому они перекомпилируются при каждой пересборке модели. Препроцессор (`.F90`) включает `#ifdef X11` в `astra_main.F90` через `$(X11_FLAG)`.

### 4.5. `Makeastra`

58 строк. Вызывается с `sub=<каталог>`, `MODOBJS=<список>`, `OBJS=<список>`, `INC_FLAGS=<-I...>`.

- Пути: `src/$(sub)` → `obj/$(sub)`, `mod/$(sub)`, `lib/lib$(sub).a`, `INC_DIR=$(AWD)/inc`.
- Если задан `SPIDER_LIB` — `SPIDER_FLAG=-DSPIDER` (строки 13–18), передаётся только для `.F90`.
- Цели: `modules` = `directories $(MODOBJS) mvmod`; `lib` = `ar rcs lib<sub>.a $(MODOBJS) $(OBJS)`; `mvmod` = удалить `*genmod.*`, перенести `*.mod` из текущего каталога (`$AWD`) в `mod/<sub>/`; `clean`.
- Правила (строки 48–55): `.f90` → `$(FC) $(FC_FLAGS) $(INC_FLAGS) -c`; `.F90` → то же плюс `$(SPIDER_INC) $(SPIDER_FLAG)`; `.c` → `$(CC) -fPIC -c ... -I inc`.

Компилятор кладёт `.mod` в текущий каталог (`$AWD`), а `mvmod` их переносит. Поэтому в корне установки во время сборки ненадолго появляются `*.mod`, а `make clean` удаляет оставшиеся (`rm -f *.mod`).

### 4.6. `Makesub`

107 строк. Компилирует **все** `*.f90`, `*.F90`, `*.f`, `*.c` из `$(AWD)/$(sub)` (`sbr` или `fnc`) в `lib/lib$(sub).a`; модули — в `mod/$(sub)`.

Флаги наличия внешних библиотек (строки 15–41) — от переменных из `astra_rc`:

| переменная | флаг препроцессора | файл в `sbr/` |
|---|---|---|
| `TORIC_LIB` | `-DTORIC_INSTALLED` | `torfpql_mod.F90` |
| `TORBEAM_LIB` | `-DTORBEAM_INSTALLED` | `torba.F90` |
| `RABBIT_LIB` | `-DRABBIT_INSTALLED` | `rabbit.F90` (+ `-I$(RABBIT_INC)`) |
| `QLKNN_LIB` | `-DQLKNN_INSTALLED` | `qlknn_serial.F90` (+ `-I$(QLKNN_INC)`) |
| `TGLF_LIB` | `$(TGLF_FLAG)` — **не определяется** | `tglf_serial.F90` (`#ifdef TGLF_INSTALLED`) |

Без флага файл компилируется в заглушку (`#else`-ветка; коммит `2ec99235`, май 2026). **Находка**: в `Makesub` нет блока `ifdef TGLF_LIB → TGLF_FLAG = -DTGLF_INSTALLED`, поэтому `tglf_serial.F90` по коду всегда собирается заглушкой, даже если TGLF установлен (если только `TGLF_FLAG` не задан в окружении). Не проверено сборкой с TGLF.

Особые правила (строки 79–95): `tglf_serial.o`, `qlknn_serial.o`, `torba.o`, `rabbit.o`, `torfpql_mod.o` (с `-I$(MOD_DIR)`), `nbi.o` (с `-I mod/nbi`). Общие правила: `-I mod/misc -I mod/astra`. Цель `all` печатает список исходников (`print`), поэтому в логе видно `echo "dir ..."`.

### 4.7. `Makexpr`

84 строки. Собирает три **отдельные программы** для параллельных расчётов через System V IPC (разделяемая память + семафоры; запускаются из `src/astra/ipc_control.c` через `system("<AWD>/<имя> ... &")`):

| цель | исходники | линковщик | библиотеки |
|---|---|---|---|
| `xpr/tglfi` | `xpr/tra_interf.c` (`-Dtglf`) + `xpr/tglf_interf.f90` | `$(FC)` (gfortran) / `$(CC)` (Intel) | `lib/libastra.a`, `$(TGLF_LIB)`, MKL |
| `xpr/qlki` | `tra_interf.c` (`-Dqlk`) + `qlk_interf.f90` (через `$(MPIFC)`) | `$(MPIFC) -fopenmp` / `$(MPICC) -qopenmp` | `libastra.a`, `$(QLK_LIB)`, MKL |
| `xpr/neo` | `tra_interf.c` (`-Dneo`) + `neo_interf.f90` (через `$(MPIFC)`) | как для `qlki` | `libastra.a`, `$(NEO_LIB)`, MKL |

- Ветка компилятора (строки 12–23): при `FC=gfortran` — OpenMP `-fopenmp`; иначе Intel-библиотеки `-lifcoremt -lifport -lifcore` из `COMP_LIBDIR`.
- `directories` создаёт рабочие каталоги `$(AWD)/tglf`, `$(AWD)/qualikiz`, `$(AWD)/neo`.
- Цели зависят от `exe/astra_rc` и `lib/libastra.a`, так что пересобираются при изменении ядра.
- Без библиотек линковка падает, `Build` выходит с кодом 2 **до запуска модели**. Upstream убирает цели только в `install.sh -safe`; наш установщик — всегда ([`DETAILED.md`](../../installer/docs/DETAILED.md), разд. 3.6). В `~/astra/a8/exe/Makexpr` строка 36 после патча: `all: directories`.

### 4.8. `astra_rc`, `astra_rc_gcc`, `astra_rc_ifx`

Задают `*_HOME`, `*_LIB`, `*_INC` внешних модулей внутри `$ASTRA_EXT`. Принцип: если библиотека-файл существует — переменная `*_LIB` экспортируется, иначе `unset` (тогда `Makesub` собирает заглушку, а `Makefile` не линкует библиотеку). Исключения: `QLK_LIB` и `NEO_LIB` экспортируются **без проверки** существования.

`astra_rc` (используется по умолчанию) — версии каталогов (на `67281112`):

| модуль | путь в `$ASTRA_EXT` (Linux) | на `darwin` | файл-признак | переменные |
|---|---|---|---|---|
| json-fortran | `json/9.0.2` | `json` | — (не проверяется!) | `JSON_LIB=.../lib/libjsonfortran.a`, `JSON_INC=.../inc` |
| STELLOPT | `STELLOPT/` (если `STELLOPT_HOME` не задан) | то же | — | `STELLOPT_HOME` |
| SPIDER | `spider/latest` | то же | `lib/libspider.a` | `SPIDER_LIB`, `SPIDER_INC=-I.../inc` |
| TORBEAM | `torbeam/nov24` | `torbeam` | `lib/libtorbeamB.so` | `TORBEAM_LIB` (+ `-lirc` из `COMP_LIBDIR`) |
| RABBIT | `rabbit/may26` | `rabbit` | `lib/librabbit.so` или `.dylib` | `RABBIT_LIB`, `RABBIT_INC` |
| TORIC | `toric/apr26/intel/2025.1` | то же | `lib/libtoric.so` | `TORIC_LIB` (FFTW по жёсткому пути `/mpcdf/soft/SLE_15/...`) |
| QuaLiKiZ | `qualikiz/nov24` | `qualikiz` | — (не проверяется) | `QLK_LIB` (intel, или gcc если есть), `QLK_INC` |
| QLKNN | `qlk_nn/nov24` | `qlk_nn` | `libQLKNN-{gcc,intel}-release-default.a` | `QLKNN_LIB`, `QLKNN_INC` |
| TGLF | `tglf/jun25` | `tglf` | `lib/libtglf.a` | `TGLF_LIB`, `TGLF_INC` |
| NEO | `neo/dec24` | `neo` | — (не проверяется) | `NEO_LIB` (6 архивов), `NEO_INC` |

`astra_rc_gcc` и `astra_rc_ifx` (берутся при `ASTRA_COMPILER=gcc|ifx`, т.е. `as_exe -c gcc`) отличаются каталогами: для gcc — `json/gcc_9.0.2`, `spider/gcc_may26`, `torbeam/gcc_apr26`, `rabbit/gcc_may26`, `toric/apr26/gfortran/2025.1`, `qualikiz/gcc_apr26`, `qlk_nn/gcc_apr26`, `tglf/gcc_apr26`, `neo/gcc_apr26`; для ifx — `json/9.0.2`, `spider/may26`, `torbeam/nov24`, `rabbit/may26`, `qualikiz/nov24`, `qlk_nn/nov24`, `tglf/jun25`, `neo/dec24`. В них нет ветки `darwin` и строки `STELLOPT_HOME`. Пары `env.<pl>_gcc`/`env.<pl>_ifx` есть только для `tok` (`platform/env.tok_gcc`, `platform/env.tok_ifx`).

Практическое следствие: **json-fortran — единственная обязательная внешняя библиотека**, и путь к ней не проверяется: при её отсутствии ошибка появится только при компиляции (`json_module` не найден) или линковке. На Ubuntu установщик кладёт её ровно в `$ASTRA_EXT/json/9.0.2` ([`DETAILED.md`](../../installer/docs/DETAILED.md), разд. 2).

### 4.9. `greenMatrices.py`

600 строк. Считает матрицы Грина и индуктивности катушек и пассивных проводников по описанию машины (эллиптические интегралы `scipy.special.ellipk/ellipe`, `mu_0`) — вход для FEQIS в режиме свободной границы.

`main()` (строки 572–596): для каждого `exp/cnf/*.json` берёт префикс до `_` (`aug`, `iter`, `jet`), и если есть `exp/cnf/<tok>_description_in.json`, пишет `exp/cnf/<tok>_description_full.json`, когда выход отсутствует, пустой, старше входа или старше самого `greenMatrices.py`. Если во входе нет `Rmin`, расчёт пропускается с предупреждением (так для `iter` и `jet` в поставке: `"Rmin" not found, skipping Green-function calculation`).

На стенде: `exp/cnf/aug_description_full.json` — 11.7 МБ, расчёт ~1 с (`~/astra/test_flux_feqis.log`, строки 1–13: `nactive, ncoils, nconduc 12 21 52`). Логгер `green` выставлен в `DEBUG`, поэтому эти строки есть в выводе каждого запуска, когда файл пересчитывается. После обновления ASTRA (новый `greenMatrices.py`) первый запуск всегда пересчитывает файл.

Описание машины читают: `src/astra/machine_config.f90::config_read` (`exp/cnf/<machine>_description_in.json`) и `src/feqis/feqis_solvers.f90:561` (`..._description_full.json`). Подробно — [06-equilibrium.md](06-equilibrium.md).

### 4.10. `version`, `wr_nml`, `README`

- **`version`** — баннер «ASTRA - Automated System of TRansport Analysis, Version 8.6, January 2026», создатели и разработчики. Читается графикой (`src/graph/graph_utils.f90:172`: `open(131, FILE='exe/version')`) для заголовка окна «ASTRA 8.6 -- Model: ... -- Data: ...». Чтобы поменять отображаемую версию, достаточно поправить этот файл.
- **`wr_nml`** — 14 строк: `source /etc/profile.d/modules.sh`, `module load trview`, `python ${TRVIEW_HOME}/trview/astra_nml_u.py "$@"`. Это внешний инструмент IPP (TRVIEW, по названию — запись namelist для ASTRA); вне IPP не работает. Назначение по коду не установить — сам `astra_nml_u.py` в репозиторий не входит.
- **`README`** — 10 строк: перечисляет `mod/{misc,feqis,fnc,astra,nbi,sbr}` и `lib/lib*.a`; говорит, что SPIDER, TGLF, QuaLiKiZ, QLKNN, json, NEO предкомпилированы и общие для площадки, а версию выбирают в `exe/astra_rc` (в тексте опечатка `$AWS/exe/astra_rc`).

## 5. Внешние модули и `$ASTRA_EXT`

### 5.1. Что ищется и где

`ASTRA_EXT` задаётся в `platform/env.<pl>` (например, `/tokp/work/software/TOK_2023/software/ASTRA_LIBRARIES_EXT` в `env.tok`, `${AWD}/astra_modules` в `env.darwin`, `$HOME/astra/ext` в нашем `env.ubuntu`). Раскладка внутри — по `exe/astra_rc` (таблица в 4.8). Программа также читает `ASTRA_EXT` во время выполнения (`src/astra/read_input.f90:68`, `getenv`), предположительно для поиска данных модулей.

Раскладка на стенде (`~/astra/ext`): только `json/9.0.2/{lib,inc}` — этого достаточно для сборки и теста.

### 5.2. `modules_install.sh` (кластеры)

515 строк, интерактивный: для каждого пакета `read -p "Install X (y/n)"`. Аргумент — необязательный суффикс компилятора (`gcc`/`ifx`): площадка становится `<pl>_<comp>`. Сначала `source platform/env.<pl>`; исходники качаются в `$HOME/soft`.

| пакет | откуда | куда в `$ASTRA_EXT` | доступ |
|---|---|---|---|
| Intel oneAPI HPC 2025.2 | registrationcenter-download.intel.com, ~2 ГБ | `$HOME/soft/intel/2025.02/oneapi` | открыто |
| CMake 3.30.3 | GitHub Kitware | `$HOME/soft/cmake-3.30.3-linux-x86_64` | открыто |
| json-fortran 9.0.2 | GitHub jacobwilliams | `json/9.0.2` (Intel) или `json/gcc_9.0.2` (gfortran) | открыто |
| NetCDF (zlib 1.3.2, HDF5 1.14.3, netcdf-c 4.9.2, netcdf-fortran 4.6.1) | zlib.net, hdfgroup, GitHub Unidata | `netcdf/mar26` или `netcdf/gcc_mar26` | открыто |
| RABBIT | `git@gitlab.mpcdf.mpg.de:markusw/rabbit.git` | `rabbit/unstable` | **IPP, SSH** |
| TORBEAM | `git@gitlab.mpcdf.mpg.de:ipp-aug/torbeam` | `torbeam/unstable` | **IPP, SSH** |
| SPIDER | `git@gitlab.mpcdf.mpg.de:git/spider` | `spider/unstable` | **IPP, SSH** |
| QuaLiKiZ | `https://gitlab.com/qualikiz-group/QuaLiKiZ.git` | `qualikiz/unstable` | открыто, нужен MPI |
| QLKNN + 4 набора namelist'ов | `https://gitlab.com/qualikiz-group/QLKNN-fortran.git` и др. | `qlk_nn/unstable` | открыто |
| TGLF, NEO (GACODE) | `git@github.com:gafusion/gacode.git` (SSH!) | `tglf/unstable`, `neo/unstable` | открыто, но клон по SSH требует ключа GitHub; нужен MPI |
| STRAHL | `git@gitlab.mpcdf.mpg.de:rld/strahl.git` | `strahl/unstable/bin` | **IPP, SSH**; `machine`-файл с `ifort` зашит |

В каждый каталог установки копируется `platform/env.<pl>` (какой средой собрано) и файл `hash` с коммитом исходников. В конце — `chmod -R a+rx $ASTRA_EXT` (общий доступ на площадке).

Несоответствия, найденные при чтении:

- Модули ставятся в каталоги `.../unstable`, а `exe/astra_rc` ищет `nov24`, `may26`, `jun25`, `dec24`, `latest`. После установки нужно либо переименовать/сделать симлинки, либо править `astra_rc`. Для gfortran json-fortran уходит в `json/gcc_9.0.2`, а `astra_rc` по умолчанию ищет `json/9.0.2` (совпадает только с `astra_rc_gcc`).
- QuaLiKiZ: всегда копируется `libQLK-intel-release-default.a`, даже при `gcc` (строка 323).
- `NPROC` используется в `make -j$NPROC` (сборка NetCDF), но нигде не задан — `make -j` без числа (неограниченный параллелизм).
- Строка 200: `export FFLAGS=FFLAGS_IN` — присваивается буквальная строка вместо `$FFLAGS_IN`.
- Шаг установки самой ASTRA выключен (`ASTRA_FLAG=n`, строка 503).

### 5.3. `modules_inst_loc.sh` (локальные машины)

1043 строки (коммит `7698c2f7` «More portability for MacOS, Linux PC», 2026-08-10). Две ветки: `platform == "linux"` (apt) и `platform == "darwin"` (Homebrew). Раскладка **плоская**: `$ASTRA_EXT/{json,netcdf,netcdf_fortran,rabbit,torbeam,qualikiz,qlk_nn,tglf,neo,strahl}` — она совпадает с веткой `darwin` в `astra_rc`. SPIDER эта версия не ставит. TORBEAM берётся из другого места: `ssh://gerrit.ipp.mpg.de:29418/ipp/e1/ecrh/libtorbeam`, STRAHL — `git@gitlab.mpcdf.mpg.de:azito/strahl.git`.

Linux-ветка: проверяет/предлагает поставить через apt `gcc`, `cmake` (обязательно), `libomp-dev`, `libopenmpi-dev`, `libopenblas-dev`; при CMake < 3.18 качает 3.30.3 (x86_64/aarch64); json-fortran собирается в `/tmp` и ставится в `$ASTRA_EXT/json`; NetCDF не собирается, а **копируется** из системных пакетов `libnetcdf-dev`, `libnetcdf-fortran-dev`/`libnetcdff-dev` в `$ASTRA_EXT/netcdf`, `netcdf_fortran`.

**Находка**: ветка Linux на практике недостижима. `get_platform` никогда не печатает `linux` (обычный ПК — `mypc`), файла `platform/env.linux` нет, а суффикс компилятора добавляется, только если есть `env.<pl>_<comp>`. Значит, на Linux скрипт либо падает на `source platform/env.mypc`, либо (с нашим `env.mypc` → `env.ubuntu`) проходит мимо обеих веток и ничего не ставит. Кроме того, даже при успехе плоская раскладка не совпадает с Linux-веткой `astra_rc` (`json/9.0.2` и т.д.). Вывод сделан по коду, запуском не проверялся.

### 5.4. Что реально нужно на обычной Ubuntu

| модуль | нужен для | статус |
|---|---|---|
| json-fortran | любая сборка (JSON-ввод/вывод, описание машины) | обязателен; ставит `install_astra8.sh` |
| BLAS/LAPACK (MKL или OpenBLAS) | FEQIS, линковка | обязателен; OpenBLAS через `MAKEFLAGS` |
| libX11 | GUI | для окна; batch без него собирается |
| RABBIT, TORBEAM, SPIDER, STRAHL, TORIC | NBI, ECRH, SPIDER-равновесие, примеси, ICRH | закрыты (IPP); без них — заглушки, модели нужно править ([07](07-heating-and-external-modules.md)) |
| TGLF, NEO, QuaLiKiZ, QLKNN | турбулентный и неоклассический транспорт | открыты, не поставлены; нужны MPI и правка `astra_rc`/`Makexpr` |
| NetCDF | только RABBIT | нет смысла без RABBIT |
| STELLOPT (VMEC2000, BOOZ_XFORM, DKES) | стеллараторные модели (`memo_stella`) | не поставлен |

## 6. `requirements.txt` и реальные зависимости

`requirements.txt` — не файл pip, а памятка (весь текст):

```
Compiler:
   intel >2018 (ifort, ifx)
   gfortran

Additional packages
   MKL
   MPI
   OpenMP flags
   NetCDF if RABBIT is needed

python
   v >= 3
   scipy v >= 0.19.1
```

`pip install -r requirements.txt` с ним не сработает. Фактические зависимости по импортам:

| что | кто требует | обязательно |
|---|---|---|
| `python3` | `as_exe`, парсер, постобработка | да |
| `numpy` | `as_exe.py`, `greenMatrices.py`, `pyparse/*` (4 модуля), `json2cdf.py` | да |
| `scipy` | `greenMatrices.py` (`scipy.special`, `scipy.constants`) — импортируется при **каждом** запуске `as_exe.py`; `json2cdf.py`, `compareRegressions.py` (`scipy.io.netcdf_file`) | да |
| `matplotlib` | `plot_compare.py`, графические скрипты в `postProcess/`, `vmec/python/` | нет |
| `netCDF4` | не нужен: CDF пишется и читается через `scipy.io.netcdf_file` (NetCDF-3) | нет |
| `imas` | `postProcess/acdf2imas.py`, `preProcess/imas2astra.py` | нет (только IMAS) |
| IDL, csh | `strahl/pstrahl` | нет (только STRAHL на IPP) |
| `make`, `ar`, `lscpu`, `nproc`, `ipcrm` | сборка, `Build`, уборка IPC | да (`util-linux`, `binutils`) |
| MPI (`mpif90`, `mpicc`) | только `Makexpr` (QLK, NEO) | нет, если цели убраны |

Версия Python явно не ограничена; на стенде — 3.12, в `ubuntu:22.04` — 3.10 ([`DETAILED.md`](../../installer/docs/DETAILED.md), разд. 5.4).

## 7. Каталоги и файлы, которые появляются при работе

Сверено с `~/astra/a8` после установки и двух расчётов (`flux_feqis`/`aug34954`, `fbe`/`AUG33040_2500`).

| путь | кто создаёт | содержимое | `make clean` удаляет |
|---|---|---|---|
| `equ/`, `exp/`, `udb/`, `fml/`, `fnc/`, `sbr/`, `xpr/`, `exe/`, `preProcess/` | `install.sh` / установщик | копии из `repo_user_area/` | нет |
| `pyparse/iondens_ass.py` | `install.sh` | см. [01](01-overview.md), разд. 6.3 | нет |
| `tmp/` | `install.sh` (`mkdir -p tmp`) | `astra.nml`; при IPC — `<e><m>-<id>.ipc` | нет |
| `exp/cnf/<dev>_description_full.json` | `greenMatrices.py` | матрицы Грина (aug — 11.7 МБ) | нет |
| `exe/__pycache__/` | Python при импорте `greenMatrices` | байткод | нет |
| `bin/` | `Makefile` (`directories`) | `<m>_gui.exe`, `<m>_batch.exe` | да |
| `src/tmp/` | `Makefile` + парсер | 9 `.f90`, `model.txt`, `declar.fml`, `declar.fnc` — для **последней** собранной модели | нет |
| `lib/` | `Makefile`, `Makeastra`, `Makesub` | `libastra.a`, `libfeqis.a`, `libfnc.a`, `libgraph.a`, `libmisc.a`, `libnbi.a`, `libsbr.a` | да |
| `mod/{astra,feqis,fnc,graph,misc,nbi,sbr}/` | `Makeastra`, `Makesub` | `.mod` | да |
| `obj/{astra,feqis,fnc,graph,misc,nbi,sbr}/` | то же | `.o` | да |
| `ncdf_out/` | `Makefile` (`directories`), программа, `json2cdf` | `<e><m>-<n>.json`, `<e><m>.CDF` | нет |
| `dat/` | `Makefile` (`directories`) | рабочие файлы модулей (FEQIS: `dat/feqis_time{1,2}_file.dat`, `dat/output_timefit.dat`, ...); на стенде пуст | нет |
| `tglf/`, `qualikiz/`, `neo/` | `Makexpr` (`directories`) | рабочие каталоги IPC-программ; на стенде пусты | нет |
| `xpr/*.o`, `xpr/tglfi`, `xpr/qlki`, `xpr/neo` | `Makexpr` | IPC-программы | да (`Makexpr clean`) |
| `exe/batch.sh` | `Build` при `sbatch` | скрипт задания SLURM | нет |
| `regressions/<e><m>.log` | `compareRegressions.py` | отчёт сравнения | нет |
| `*_backup/` | `install.sh -safe` | резервные копии каталогов | нет |
| `vmec/dat` и др. | стеллараторная часть | рабочие файлы VMEC (игнорируются git) | нет |

## 8. Полная последовательность одного запуска

Пример: `astra8 -m flux_feqis -v aug34954 -s 4 -e 5 -batch` на стенде (24 ядра, WSL2). Номера строк лога — `~/astra/test_flux_feqis.log` (369 строк).

1. **`astra8`** (`~/.local/bin/astra8`, создаёт установщик): `export ASTRA_PLATFORM=ubuntu`, `export AWD=~/astra/a8`, `cd $AWD`, `exec exe/as_exe "$@"`.
2. **`exe/as_exe`**: `get_platform` → `ubuntu`; `source platform/env.ubuntu` (`FC=gfortran`, `FC_FLAGS`, `LD_FLAGS`, `MAKEFLAGS`, `OPENBLAS_NUM_THREADS=1`, `ASTRA_EXT`); `PATH` с `/usr/bin`; `python3 exe/as_exe.py ...`.
3. **`greenMatrices`** (лог 1–13): пересчёт `aug_description_full.json` при необходимости; пропуск `jet`, `iter`.
4. **`tmp/astra.nml`** (лог 14): `Writing /home/egor/astra/a8/tmp/astra.nml`.
5. **`exe/Build flux_feqis aug34954 BGD False 1 .../stepup.F90`** (лог 15–16). Окружение ещё раз: `env.ubuntu`, `exe/astra_rc`; `ulimit -s`; `X11` снят, `ATASK=bin/flux_feqis_batch.exe`.
6. **`make -f exe/Makefile`** (лог 17–310):
   - `modules` — 5 вызовов `Makeastra ... modules` (misc, astra, feqis, graph, nbi);
   - `lib*.a` — `Makeastra ... lib` для каждого `src/<sub>`;
   - `Makesub sub=sbr` (лог 117–164), `Makesub sub=fnc` (165–292);
   - парсер: `pyparse/parser_main.py -equ equ/flux_feqis -exp exp/aug34954` (лог 293–308) — `Written file ./src/tmp/...`, `NRX_MAX 51`;
   - линковка (лог ~309): `gfortran -fopenmp -Wl,--no-as-needed -I$AWD -o bin/flux_feqis_batch.exe src/astra/astra_main.F90 src/astra/stepup.F90 src/tmp/{associate_pointers,inivar,detvar,postep,ininam,init_converge_step,eqns_inc}.f90 -I mod/... -lsbr -lnbi -lastra -lfnc -lfeqis -lmisc ext/json/9.0.2/lib/libjsonfortran.a -L lib -lopenblas -lrt`. Предупреждение компоновщика `json_value_module.F90.o: requires executable stack` безвредно.
7. **`make -f exe/Makexpr`** (лог 311): `Compile IPC programs`; цели убраны — ничего не делает.
8. **Проверка конфликтов** `libsbr.a`/`libfnc.a` — молча.
9. **Запуск** (лог 312): `Running .../bin/flux_feqis_batch.exe without graphic interactive window`. Программа: `Note: The following floating-point exceptions are signalling: IEEE_INVALID_FLAG ...` (сообщение gfortran при `STOP`, не ошибка), чтение U-файлов `udb/34954/*`, границы, итерации FEQIS (`surf2surf time`), `Written file .../ncdf_out/aug34954flux_feqis-<n>.json` каждые `DPOUT`.
10. **Выход** (лог 353–359): `>>> ASTRA normal exit >>>`, `Total time steps 21`, `Average time step 50.000 msec`, `Astra wall time 00:00:30`, разбивка по подпрограммам (`Equilibrium 98.8%`).
11. **Уборка IPC**: файлов `tmp/*.ipc` нет — ничего.
12. **`json2cdf`** (лог 361–369): 6 JSON → `ncdf_out/aug34954flux_feqis.CDF` (`nt=6, nrho=91, nr_eq=91, nthe_eq=90, nR=65, nZ=65`).

Время по замерам ([`DETAILED.md`](../../installer/docs/DETAILED.md), разд. 5.1): полная сборка с `-O0` ~10 с (205 файлов; make вызывается без `-j`), расчёт 30 с. Повторный запуск той же модели без изменений: `make` ничего не пересобирает, сразу запуск.

`json2cdf` склеивает файлы `-1`, `-2`, ... по порядку и останавливается на первом отсутствующем номере или на файле, который **старше** предыдущего (`postProcess/json2cdf.py::json_concat`, строки 42–51). Так отбрасываются хвосты от прежнего, более длинного расчёта с тем же `<e><m>`. Если `-1.json` нет, функция просто возвращается, а `as_exe.py` подавляет любые исключения — CDF не появится, и сообщения об этом не будет.

## 9. Batch, GUI, SLURM, flight simulator

| режим | как включить | исполняемый файл | `TASK` | что происходит |
|---|---|---|---|---|
| GUI (активный) | без `-batch` | `bin/<m>_gui.exe` (`-DX11`, `libgraph.a`, libX11, `astra_out.f90`, `set_graph_names.f90`) | `RUN` | окно X11 «ASTRA 8.6 -- Model: ... -- Data: ...»; пауза в `-tpause`; программа **не завершается сама** после `TEND` — ждёт действий пользователя (по наблюдениям, [`../../AGENTS.md`](../../AGENTS.md)) |
| batch локально | `-batch`, `sbatch` нет в `PATH` | `bin/<m>_batch.exe` (без X11) | `BGD` | расчёт до `TEND`, выход |
| batch в SLURM | `-batch`, есть `sbatch` | то же | `BGD` | `exe/batch.sh` + `sbatch [-W]`; `as_exe.py` после `sbatch` сразу пытается склеить JSON (без `-W` — ещё не готовых) |
| flight simulator | `-fs` | `bin/<m>_batch.exe` | `RUN`/`BGD` | `flightsim=.True.`; режимы `IPEQL=-1/-6` обходятся, данные должны приходить от внешней системы управления через подпрограммы в `sbr/` (`wrsudg.f90`, `rrsudg.f90`, `shmrw.c` — в репозиторий не входят; `docu/addendum_A8.tex`, разд. 5) |

GUI-сборка и batch-сборка одной модели — **разные** исполняемые файлы с общими библиотеками, поэтому переключение режимов вызывает линковку (и парсер), но не пересборку библиотек. Окно на WSL2/Windows 11 выводит WSLg; без X-сервера — `-batch` ([`DETAILED.md`](../../installer/docs/DETAILED.md), разд. 5.2 и 8). Устройство окна — [08-graphics-and-output.md](08-graphics-and-output.md).

## 10. Пересборка и очистка

### 10.1. Что вызывает пересборку

По зависимостям в `exe/Makefile` (строка 210 и цели библиотек):

| изменили | пересобирается | парсер перезапускается |
|---|---|---|
| `equ/<m>` | исполняемый файл модели | да |
| `fml/*`, `pyparse/*.py` | исполняемый файл **текущей** модели | да |
| `exe/astra_rc` | исполняемый файл; также IPC-программы | да |
| `sbr/*` | `libsbr.a` (только изменённые `.o`), затем исполняемый файл | да |
| `fnc/*.f90` | `libfnc.a`, исполняемый файл | да |
| `src/<sub>/*` | `lib<sub>.a` (изменённые объекты), исполняемый файл | да |
| `src/astra/astra_main.F90`, файл `-workflow` | исполняемый файл | да |
| `exp/<e>` | **ничего** | нет |
| `equ/log/<m>`, `exp/nml/*`, `udb/*`, `exp/cnf/*` | ничего (читаются во время выполнения) | нет |
| `platform/env.*` (компилятор, флаги) | **ничего** | нет |
| переход на другую модель | только исполняемый файл новой модели (если его нет или он устарел) | да |

Следствия:

- **Смена флагов компилятора или BLAS** не видна make: после правки `env.<pl>` нужен `./clean.sh`.
- **Зависимости модулей Fortran не отслеживаются**: при изменении интерфейса модуля в `src/` перекомпилируется только изменённый файл, а зависящие от него — нет. Надёжнее `./clean.sh`.
- **Файл эксперимента не является зависимостью**, хотя из него берётся `nr_x_max` — размер массивов `XAXES`, `DATAX` (`pyparse/parser_main.py:30`, `pyparse/code_gen.py:98`, `src/astra/read_input.f90:62`). Если запустить уже собранную модель с другим экспериментом, где у профилей больше радиальных точек, возможен выход за границы массива. Вывод сделан по коду и запуском не проверялся; безопасный путь — `touch equ/<m>` перед расчётом с новым экспериментом.
- **Путь установки** вшит в исполняемый файл (`src/tmp/ininam.f90`: `AWD = "<путь>"`): после переноса или копирования каталога установки надо пересобрать.

### 10.2. Как чистить

| задача | команда |
|---|---|
| полная пересборка | `./clean.sh` (= `make -f exe/Makefile clean`: `*.mod` в корне, `xpr/*.o` и IPC-программы, `mod/`, `bin/`, `obj/`, `lib/`); следующий `as_exe` соберёт всё заново |
| пересобрать одну модель | `touch equ/<m>` или `rm bin/<m>_*.exe` |
| только библиотеку `sbr` | `rm lib/libsbr.a` (или `touch` изменённого файла) |
| сбросить «последний расчёт» | `rm tmp/astra.nml` (тогда `-m`, `-v` снова обязательны) |
| убрать результаты | `rm ncdf_out/<e><m>*` — особенно если новый расчёт короче старого (см. раздел 8 про `json2cdf`) |
| пересчитать матрицы Грина | `rm exp/cnf/<dev>_description_full.json` |
| обновить `exe/`, примеры из репозитория | повторить установщик (`install_astra8.sh --no-test`) или вручную `cp` из `repo_user_area/` (`install.sh` перезапишет и ваши правки) |

`clean.sh` не трогает `src/tmp/`, `tmp/`, `ncdf_out/`, `dat/` и пользовательскую область.

### 10.3. Сборка без запуска (вручную)

По `exe/Build`, строки 16–49 — так можно собрать модель, не запуская её (команда составлена по коду и в этой главе не выполнялась):

```bash
cd ~/astra/a8
export AWD=$PWD ASTRA_PLATFORM=ubuntu
source platform/env.$(./get_platform)
source exe/astra_rc
make -f exe/Makefile MODEL=flux_feqis FEXP=aug34954 \
     ATASK=$AWD/bin/flux_feqis_batch.exe WORKFLOW=$AWD/src/astra/stepup.F90
# для GUI: export X11=y и ATASK=$AWD/bin/flux_feqis_gui.exe
```

## 11. Переменные окружения

| переменная | кто задаёт | кто читает | смысл |
|---|---|---|---|
| `AWD` | `install.sh`, `Build` (export), `astra8`; в `as_exe` — локально | все `Make*` (обязательна), `wr_nml`, `astra_rc` | корень установки |
| `ASTRA_PLATFORM` | `astra8` (наш патч) | `get_platform` (только после патча) | принудительное имя площадки |
| `ASTRA_COMPILER` | `as_exe.py -c`, пользователь | `get_platform` (суффикс), `Build` (`astra_rc_<comp>`) | вариант компилятора |
| `ASTRA_EXT` | `env.<pl>` | `astra_rc*`, `modules_*.sh`, программа (`read_input.f90:68`), `batch.sh` | корень внешних модулей |
| `FC`, `CC`, `MPIFC`, `MPICC` | `env.<pl>` | `Make*` | компиляторы |
| `FC_FLAGS`, `LD_FLAGS` | `env.<pl>` (+ `-I$AWD` в `Makefile`) | `Make*` | флаги |
| `MKL_LIBDIR`, `COMP_LIBDIR`, `OPENBLAS_LIBDIR` | `env.<pl>` | `Makefile`, `Makexpr`, `astra_rc` (`-lirc`) | пути к библиотекам |
| `MAKEFLAGS` | `env.ubuntu` (`MKL_LIB=-lopenblas`) | GNU make | перекрытие `MKL_LIB` |
| `OPENBLAS_NUM_THREADS` | `env.ubuntu` | OpenBLAS | 1 поток BLAS внутри OpenMP |
| `PYTHON_BINDIR`, `PYTHON_BIN` | `env.<pl>` | `as_exe`, `compReg.sh`; `a2vmec.f90:41` (`PYTHON_BIN`) | Python для скриптов и стеллараторной части |
| `X11` | `Build` | `Makefile` | GUI-сборка |
| `N_NODES` | `Build` (из `-nodes`) | `a2vmec.f90:313` | число узлов для параллельных задач VMEC/DKES |
| `MAX_NWORKERS` | `Build` (GUI: `nproc`/потоки на ядро; SLURM: 128; локальный batch: не задаётся) | `sbr/*_ipc.f90` | число IPC-процессов |
| `SL_QUEUE` | `env.<pl>` | `Build` → `batch.sh` | директивы `#SBATCH` |
| `STELLOPT_HOME`, `MPI_COMMAND` | `env.tok`, `astra_rc` (значение по умолчанию) | `a2vmec.f90:40–42` | стеллараторная цепочка |
| `JSON_*`, `SPIDER_*`, `TORBEAM_LIB`, `RABBIT_*`, `TORIC_*`, `QLK_*`, `QLKNN_*`, `TGLF_*`, `NEO_*` | `astra_rc*` | `Makefile`, `Makeastra`, `Makesub`, `Makexpr` | внешние модули |
| `OMP_NUM_THREADS` | только в `SL_QUEUE` площадок (например, `env.tok`) | OpenMP | локально не задаётся — OpenMP берёт все логические ядра |
| `LD_LIBRARY_PATH` | `env.<pl>` (Intel, darwin) | динамический загрузчик, `batch.sh` | `.so` внешних модулей (у RABBIT/TORBEAM ещё и `rpath`) |
| `TRVIEW_HOME` | `module load trview` | `wr_nml` | только IPP |

## 12. Типичные ошибки

Проблемы, характерные для Ubuntu (CRLF, нет `env.mypc`, длина строк gfortran, `XOpenDisplay`, MKL, `Makexpr`, закрытые модули в тестовой модели, U-строки A7, код возврата `as_exe`), подробно разобраны в [`../../installer/docs/DETAILED.md`](../../installer/docs/DETAILED.md), разделы 3 и 8. Ниже — то, что вытекает из устройства конвейера.

| симптом | причина | что делать | источник |
|---|---|---|---|
| `as_exe.py: error: the following arguments are required: -m/--equ, -v/--exp` | первый запуск, нет `tmp/astra.nml` | указать `-m` и `-v` | `as_exe.py:46–48` |
| `$(error setenv/export AWD first)` | ручной `make` без `AWD` | `export AWD=$PWD` | `Makefile:1–3` |
| `platform/env.<pl>_gcc: No such file` | `-c gcc`/`ifx` на площадке без пары `env.<pl>_<comp>` (есть только для `tok`) | не использовать `-c` или создать файл | `get_platform:72–75` |
| `Error from exe/Makefile-compilation, exiting` | ошибка компиляции/линковки/парсера | искать первую `Error` выше; код возврата всё равно 0 | `Build:45–49` |
| `Error from exe/Makexpr-compilation, exiting` | нет TGLF/QLK/NEO или MPI | убрать цели из `exe/Makexpr` (строка 36) | `Build:53–58`; [`DETAILED.md`](../../installer/docs/DETAILED.md), 3.6 |
| `Conflicting X.o: in libsbr.a and libfnc.a`, `No ASTRA simulation started` | одноимённые файлы в `sbr/` и `fnc/` | переименовать один из них | `Build:62–79` |
| `Fatal Error: Cannot open module file 'json_module.mod'` или `undefined reference to json_...` (по смыслу; точный текст не воспроизводился) | нет json-fortran по пути `exe/astra_rc` | проверить `$ASTRA_EXT/json/9.0.2/{lib,inc}` | `astra_rc:7,23–24` |
| `>>> Error: file "equ/log/<m>" missing` | у модели нет файла параметров (например, `showdata`) | создать `equ/log/<m>` (можно скопировать от похожей модели) | `read_input.f90:99–103` |
| `>>> read_equ_log: Error, empty model file name` | пустое `equ_file` в `tmp/astra.nml` | указать `-m` | `read_input.f90:86–89` |
| `FileNotFoundError: .../udb/...` (парсер) или ошибка чтения U-файла (программа) | файл эксперимента ссылается на отсутствующий U-файл | положить файл в `udb/` | [`DETAILED.md`](../../installer/docs/DETAILED.md), разд. 8 |
| `X array <V> used in equ, but missing in exp` (предупреждение парсера) | модель использует `<V>X`, которого нет в эксперименте | добавить данные или изменить модель | `parser_main.py:20–22` |
| после правки флагов в `env.<pl>` ничего не изменилось | `make` не видит изменений окружения | `./clean.sh` | раздел 10 |
| странные ошибки после правки модуля в `src/` | не пересобраны зависящие модули | `./clean.sh` | раздел 10 |
| расчёт с новым экспериментом падает или даёт мусор, а со старым работал | `nr_x_max` от старого эксперимента (предположительно) | `touch equ/<m>` | раздел 10 |
| исполняемый файл пишет JSON в старый каталог / не находит файлы после переноса установки | путь вшит в `src/tmp/ininam.f90` | пересобрать (`touch equ/<m>` или `./clean.sh`) | раздел 10 |
| `-batch` «ничего не посчитал», вывод — `sbatch ...` | на машине есть SLURM | считать на узле очереди или убрать `sbatch` из `PATH` | `Build:85–107` |
| нет CDF, но и ошибки нет | исключение в `json_concat` подавлено, или JSON не записаны | `python3 postProcess/json2cdf.py -e <e><m>` вручную — увидеть ошибку | `as_exe.py:148–152` |
| CDF содержит лишние моменты времени из прошлого расчёта | `json2cdf` подхватил старые `-<n>.json` новее предыдущего | удалить `ncdf_out/<e><m>*` перед расчётом | `json2cdf.py:42–51` |
| `Segmentation fault` на больших сетках | стек | `Build` уже делает `ulimit -s 100000`; при ручном запуске сделать так же | `Build:31` |
| окно не закрывается после `TEND` | так устроен GUI-режим | закрыть окно или считать с `-batch` | [`../../AGENTS.md`](../../AGENTS.md) |
| `exe/as_exe: Permission denied` | в git `exe/*` без бита исполнения; копия сделана не через `install.sh` | `chmod 744 exe/as_exe exe/Build` | `install.sh:71–76` |

## 13. Где смотреть в исходниках

| файл | что там |
|---|---|
| `install.sh` | копирование пользовательской области, режим `-safe`, `chmod`, `make clean`, тестовый расчёт |
| `get_platform` | правила выбора площадки, суффикс `ASTRA_COMPILER` |
| `platform/env.*` | компиляторы, флаги, `ASTRA_EXT`, SLURM (см. [`platforms.md`](../../docs/astra/platforms.md)) |
| `repo_user_area/exe/as_exe` | вход: площадка, `PATH` |
| `repo_user_area/exe/as_exe.py` | опции CLI, `tmp/astra.nml`, вызов `Build`, уборка IPC, `json_concat` |
| `repo_user_area/exe/Build` | окружение, `make`, конфликты `sbr`/`fnc`, локальный запуск или `sbatch`, `MAX_NWORKERS` |
| `repo_user_area/exe/Makefile` | список генерируемых файлов, X11, порядок модулей, линковка, зависимости, `clean` |
| `repo_user_area/exe/Makeastra` | правила для `src/<sub>`, `-DSPIDER`, перенос `.mod` |
| `repo_user_area/exe/Makesub` | сборка `sbr/`, `fnc/`, флаги `*_INSTALLED` |
| `repo_user_area/exe/Makexpr` | IPC-программы `tglfi`, `qlki`, `neo` |
| `repo_user_area/exe/astra_rc`, `astra_rc_gcc`, `astra_rc_ifx` | раскладка `$ASTRA_EXT`, проверка наличия библиотек |
| `repo_user_area/exe/greenMatrices.py` | матрицы Грина для FEQIS, правило пересчёта |
| `repo_user_area/exe/version`, `wr_nml`, `README` | баннер версии, TRVIEW, справка о каталогах |
| `modules_install.sh`, `modules_inst_loc.sh` | сборка внешних модулей (кластер / локально) |
| `requirements.txt` | памятка о требованиях |
| `clean.sh` | полная очистка сборки |
| `pyparse/parser_main.py` | вход парсера, `nr_x_max`, список генерируемых файлов |
| `src/astra/astra_main.F90` | порядок работы программы, `>>> ASTRA normal exit >>>` |
| `src/astra/read_input.f90` | `read_nml` (`tmp/astra.nml`), чтение `equ/log`, `exp/nml`, `exp/<e>`, `getenv('ASTRA_EXT')` |
| `src/astra/json_rw.f90` | запись `ncdf_out/<e><m>-<n>.json`, restart |
| `src/astra/ipc_control.c`, `inc/Astra.h` | запуск IPC-программ, файлы `tmp/*.ipc` |
| `postProcess/json2cdf.py` | склейка JSON в CDF, правило остановки |
| `~/astra/a8/src/tmp/ininam.f90` | пример сгенерированного кода с вшитым `AWD` |
| `~/astra/test_flux_feqis.log` | полный лог сборки и тестового расчёта на стенде |
