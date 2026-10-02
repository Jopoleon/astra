# ASTRA 8 на Ubuntu: подробное описание установки

Документ объясняет, как устроена сборка ASTRA 8, почему штатный `install.sh` не работает на обычном Linux-ПК, что именно делает `install_astra8.sh` и как проверить результат. Короткая инструкция для пользователя — в [`../README.md`](../README.md).

Состояние на 2026-10-01:

- ASTRA upstream: `https://gitlab.mpcdf.mpg.de/git/astra`, коммит `e66f6ba1` (2026-10-01), версия в окне программы — `ASTRA 8.6`;
- стенд: Ubuntu 24.04.3 LTS в WSL2 (Windows 11, WSLg), gfortran 13.3, Python 3.12, 24 ядра;
- проверка на чистых системах: контейнеры `ubuntu:24.04` (gfortran 13.3) и `ubuntu:22.04` (gfortran 11.4), Docker 28.5;
- свежий upstream `main` (`6728111`, 2026-10-01) тоже собирается и проходит тест.

## 1. Как устроена ASTRA 8

ASTRA — транспортный код (Fortran + C + Python). В отличие от обычных программ, **исполняемый файл собирается заново под каждую модель**. Модель (`equ/<name>`) написана на собственном языке ASTRA; Python-парсер (`pyparse/`) превращает её и файл эксперимента (`exp/<name>`) в Fortran-код, который компилируется и линкуется с библиотеками ASTRA.

Цепочка вызовов:

```
exe/as_exe -m <model> -v <exp> [-s t0 -e t1] [-batch]
 └─ exe/as_exe.py                 пишет tmp/astra_log, вызывает Build
     └─ exe/Build
         ├─ get_platform           имя платформы по hostname/пути
         ├─ source platform/env.<platform>   компилятор, флаги, MKL, ASTRA_EXT
         ├─ source exe/astra_rc    пути к внешним модулям в $ASTRA_EXT
         ├─ make -f exe/Makefile   lib{misc,astra,feqis,graph,nbi,sbr,fnc}.a
         │    └─ pyparse/parser_main.py  модель+эксперимент -> src/tmp/*.f90
         │    └─ линковка bin/<model>_gui.exe или bin/<model>_batch.exe
         ├─ make -f exe/Makexpr    IPC-программы TGLF/QuaLiKiZ/NEO (MPI)
         └─ запуск bin/<model>_*.exe
 └─ json2cdf                      ncdf_out/*.json -> ncdf_out/<exp><model>.CDF
```

Первая сборка компилирует все библиотеки (2–5 минут). Повторные запуски пересобирают только сгенерированную часть модели.

Каталоги рабочей копии (`AWD`, по умолчанию `~/astra/a8`):

| каталог | назначение | git |
|---|---|---|
| `src/` | исходники ядра ASTRA | отслеживается |
| `repo_user_area/` | шаблоны пользовательской области | отслеживается |
| `equ/`, `exp/`, `udb/`, `exe/`, `sbr/`, `fnc/`, `fml/`, `xpr/`, `preProcess/` | пользовательская область, копируется из `repo_user_area/` | игнорируется |
| `platform/env.*` | настройки площадок (IPP, ITER, GA, MIT, ...) | отслеживается |
| `obj/`, `mod/`, `lib/`, `bin/` | результаты сборки | игнорируется |
| `ncdf_out/` | результаты расчёта (JSON по шагам + итоговый CDF) | игнорируется |

## 2. Внешние зависимости

| компонент | роль | источник | статус в установщике |
|---|---|---|---|
| gfortran, gcc, make, cmake | компиляция | apt | ставится |
| python3, numpy, scipy (+netCDF4, matplotlib) | парсер моделей, json2cdf, сравнение регрессий | apt | ставится |
| libx11-dev | окно с графиками | apt | ставится |
| OpenBLAS **или** MKL | LAPACK для FEQIS (свободная граница) | OpenBLAS — apt; MKL — wheel `mkl-static` с PyPI | OpenBLAS по умолчанию, `--blas mkl` — MKL |
| json-fortran 9.0.2 | ввод-вывод JSON, **обязательна** при линковке | GitHub, сборка cmake | собирается в `$PREFIX/ext/json/9.0.2` |
| NetCDF C/Fortran | нужен RABBIT | apt | ставится на будущее |
| RABBIT (NBI), TORBEAM (ECRH), SPIDER (свободная граница), STRAHL (примеси) | физические модули | **закрытые** репозитории IPP (`gitlab.mpcdf.mpg.de`, gerrit IPP) | не ставятся |
| TGLF, NEO (GACODE) | турбулентный и неоклассический транспорт | публичный GitHub `gafusion/gacode` (Apache-2.0), сборка через `mpif90` | по ключу `--with-tglf` / `--with-neo` (раздел 3.10) |
| OpenMPI (`openmpi-bin`, `libopenmpi-dev`) | компилятор `mpif90` для GACODE, линковка `xpr/neo` | apt | только с `--with-tglf`/`--with-neo` |
| QuaLiKiZ, QLKNN | турбулентный транспорт | публичные (gitlab.com), нужен MPI | пока не ставятся |

С мая 2026 (коммит `2ec99235`) ASTRA собирается без закрытых модулей: в `sbr/rabbit.F90`, `sbr/torba.F90` и в `#ifdef SPIDER` стоят заглушки. Но модели, вызывающие эти модули, нужно править: строки `TORBA...` и `RABBIT(...)` закомментировать через `!`.

## 3. Почему штатная установка не работает на Ubuntu

Найдено при сборке на стенде. Номер — как в установщике.

### 3.1. CRLF-окончания строк (WSL / Windows)

Если в git стоит `core.autocrlf=true` (так часто бывает у тех, кто работает и в Windows, и в WSL), все bash-скрипты выкачиваются с `\r\n` и не запускаются:

```
./get_platform: cannot execute: required file not found
```

Решение: `git -c core.autocrlf=false clone ...`. Установщик дополнительно ставит `core.autocrlf=false` в самом клоне.

### 3.2. Нет платформы для обычного ПК

`get_platform` угадывает площадку по hostname/пути, иначе возвращает `mypc`. Файла `platform/env.mypc` в репозитории нет, поэтому `source` падает и переменные `FC`, `MKL_LIBDIR`, `ASTRA_EXT` остаются пустыми. Второй риск — ложное срабатывание: хост с именем `node*` распознаётся как MIT, `tok*` — как IPP, `omega*` — как GA, и подтягиваются чужие пути.

Решение:

- `platform/env.ubuntu` — настройки для gfortran и библиотек из apt;
- `platform/env.mypc -> env.ubuntu` — чтобы работал и прямой вызов `exe/as_exe`;
- в начало `get_platform` добавлена переменная-переключатель `ASTRA_PLATFORM`. Команда `astra8` выставляет `ASTRA_PLATFORM=ubuntu`, так что имя компьютера больше ни на что не влияет.

Итоговый `platform/env.ubuntu` (вариант по умолчанию, OpenBLAS):

```bash
export ASTRA_EXT=$HOME/astra/ext
export FC=gfortran
export CC=gcc
export MPIFC=mpif90
export MPICC=mpicc
export FC_FLAGS="-O0 -w -fPIC -fopenmp -ffree-line-length-none"
export LD_FLAGS="-fopenmp -Wl,--no-as-needed"
export MAKEFLAGS="MKL_LIB=-lopenblas"
export OPENBLAS_NUM_THREADS=1
export PYTHON_BINDIR=/usr/bin
```

### 3.3. Строки длиннее 132 символов

```
src/astra/a2vmec.f90:251:132: Error: Invalid form of array reference at (1)
```

gfortran по умолчанию обрезает строки свободного формата на 132-м символе, ifort — нет. Решение: `-ffree-line-length-none` в `FC_FLAGS`. (На кластере IPP с gcc 14 это, по-видимому, не ловится, потому что основная сборка там идёт через Intel.)

### 3.4. Окно не линкуется: `undefined reference to XOpenDisplay`

В Ubuntu компоновщик по умолчанию работает с `--as-needed`, а в `exe/Makefile` `-lX11` стоит в строке линковки **раньше** `-lgraph` (статическая библиотека, которой X11 и нужен). В итоге libX11 выбрасывается как «ненужная». Batch-сборка при этом проходит, ломается только GUI. На SUSE (кластеры IPP) `--as-needed` по умолчанию не включён, поэтому у разработчиков ошибки нет.

Решение без правки Makefile: `-Wl,--no-as-needed` в `LD_FLAGS`.

### 3.5. MKL прописан жёстко

`exe/Makefile` и `exe/Makexpr` всегда линкуют статический MKL (`$(MKL_LIBDIR)/libmkl_intel_lp64.a ...`); OpenBLAS используется только при `PLATFORM=darwin`. Варианты:

- `--blas openblas` (по умолчанию): в `env.ubuntu` пишется `MAKEFLAGS="MKL_LIB=-lopenblas"`. Переменные в `MAKEFLAGS` GNU make воспринимает как параметры командной строки, поэтому они перекрывают присваивание внутри Makefile, в том числе во вложенных вызовах make. Сами Makefile не правятся. Дополнительно `OPENBLAS_NUM_THREADS=1`: ASTRA сама распараллелена через OpenMP, а upstream линкует **последовательный** MKL (`libmkl_sequential`). Многопоточный OpenBLAS внутри OpenMP-регионов давал переподписку: на тесте 21 мин системного времени вместо 23 с, а расчёт шёл 57 с вместо 30 с;
- `--blas mkl`: три статические библиотеки (`libmkl_intel_lp64.a`, `libmkl_core.a`, `libmkl_sequential.a`) берутся из официального wheel Intel `mkl-static` 2026.1.0 с PyPI (~250 МБ, скачался за ~1.5 мин) и кладутся в `$PREFIX/ext/mkl/2026.1.0/lib`. apt-пакет `libmkl-dev` (MKL 2020.4) тоже подходит, но весит ~1 ГБ, а зеркало Ubuntu на стенде отдавало его со скоростью ~50 КБ/с.

Выбор библиотеки на результат не влияет (раздел 5.3), поэтому по умолчанию стоит лёгкий OpenBLAS.

### 3.6. IPC-программы TGLF / QuaLiKiZ / NEO

`exe/Makexpr` по умолчанию собирает `xpr/tglfi`, `xpr/qlki`, `xpr/neo`. Для них нужны библиотеки из `$ASTRA_EXT` и MPI; `QLK_LIB` в `astra_rc` прописан без проверки существования файла. Без них `Build` завершается с кодом 2 ещё до запуска модели. Upstream-`install.sh` убирает эти цели только в режиме `-safe` и после интерактивных вопросов; установщик убирает их всегда (как ответ `n` в `-safe`). С `--with-tglf`/`--with-neo` цели `xpr/tglfi`/`xpr/neo` возвращаются (раздел 3.10).

### 3.7. Тестовая модель требует закрытых кодов

`install.sh` сразу запускает `flux_feqis`/`aug34954`, а эта модель вызывает `TORBA` (ECRH) и `RABBIT` (NBI). Установщик комментирует эти две строки в `equ/flux_feqis` (как `install.sh -safe` при ответе `n`).

### 3.8. Совместимость с ASTRA 7: файлы эксперимента и модели

Файлы и модели ASTRA 7 ломают ASTRA 8 в нескольких местах. Часть ошибок — от разницы компиляторов: A7 собиралась ifort, который прощает больше, чем gfortran. Часть — ошибки самой A8. Установщик накладывает один патч на пять файлов:

| файл | что было | что стало |
|---|---|---|
| `pyparse/exp_parser.py` | падал на U-строках с табуляцией, одномерных U-файлах, `2.0d0`, значениях за 22-й колонкой, `FACTOR 1.0d3`; тихо обрезал блоки с текстом | читает эти формы; если блок не разобран, его размер всё равно учитывается в `nr_x_max` |
| `src/astra/read_input.f90` | слово `points` в комментарии после значения открывало 2D-группу; `0.45d`, вырезанное из `0.45d0` по колонкам 17–22, gfortran не читает; для строк, не являющихся U-файлами, к размеру прибавлялся мусор | комментарии `!` обрезаются; поле расширяется до границы слова; размер считается только для U-строк |
| `src/astra/parse_utils.f90` | строка данных A7 длиннее 20 чисел → запись за конец массива, segfault | лишние слова не пишутся |
| `pyparse/parse_as.py` | `((A)*2)/(B)` → `ValueError` в генераторе (ошибка A8) | исправлено |
| `pyparse/const_text.py` | `"cc_nc(j)"` в начальных значениях не компилировался | `use nclass_mod` в `inivar` |

Проверено: вывод парсера для `aug34954`/`AUG33040_2500` и сгенерированный Fortran для всех 28 моделей A8 не изменились. Перенос моделей с формулами и `equ/log` — команда `astra8-a7to8`. Подробно — `docs/research/a7-user-area-to-a8-2026-10-02.md`.

### 3.9. `as_exe` всегда завершается с кодом 0

Даже при ошибке компиляции `exe/as_exe` возвращает 0 (результат `os.system` игнорируется). Об успехе можно судить только по выводу: `>>> ASTRA normal exit >>>` и наличию `ncdf_out/<exp><model>.CDF`. Установщик проверяет именно это.

### 3.10. TGLF и NEO из GACODE (`--with-tglf`, `--with-neo`)

TGLF (турбулентный перенос) и NEO (неоклассика) — части открытого пакета GACODE (General Atomics, Apache-2.0, <https://github.com/gafusion/gacode>). ASTRA вызывает их двумя способами:

- **IPC** (основной): обёртка `sbr/tglf_ipc.f90` (`TGLF_IPC`, `tglf_ipc(...)` в моделях `flux_tglf`, `tglf`, `tglf_pid`) через разделяемую память System V запускает 64 процесса `xpr/tglfi` — по одному на точку 64-точечной радиальной сетки — и собирает результат в `tglf_out%chi_i`, `chi_e`, `e_pflux`, ...; для NEO то же с `sbr/neo_ipc.f90` и `xpr/neo` (`flux_neo`);
- **последовательный**: `sbr/tglf_serial.F90` (`tglf_serial(...)` в `flux_tglf_serial`, `tglf_serial`) линкуется прямо в `bin/<модель>_*.exe` и считает те же 64 точки по очереди.

Что делает установщик (блок `# >>> gacode` в `install_astra8.sh`, шаг 5c):

1. apt: `openmpi-bin libopenmpi-dev` (OpenMPI 4.1.6 в 24.04, 4.1.2 в 22.04). MPI нужен только для компиляции: `tglf_init_mpi.f90` делает `use mpi`, и библиотека GACODE собирается через `mpif90`. `xpr/tglfi` линкуется обычным `gfortran` без MPI: из архива берутся только объекты без MPI-вызовов. `xpr/neo` линкуется через `mpif90` (`LD_MPI` в `Makexpr`);
2. GACODE скачивается по HTTPS одним коммитом (`git fetch --depth 1 --filter=blob:none`) и только нужные каталоги (`git sparse-checkout`: `tglf neo shared f2py/expro f2py/geo platform/build`) — ~13 МБ за 1–2 с вместо ~630 МБ полного дерева (почти всё — `cgyro`). Проверенный коммит — `master` от 2026-09-30, `f34a9dd2846afedba7a829d70f6194c2d8d8b304`; другой задаётся `--gacode-ref` (ветка, тег или полный хэш). Теги `stable_r*` в репозитории старые (последний — 2017 г.), поэтому закреплён коммит `master`. Каталог — `$PREFIX/build/gacode`;
3. платформа GACODE `platform/build/make.inc.ASTRA_UBUNTU` пишется установщиком: `FC = mpif90 -fPIC -I/-J $GACODE_ROOT/modules`, `FMATH = -fdefault-real-8 -fdefault-double-8` (как на всех площадках GACODE), `FOPT = -O2 -fallow-argument-mismatch` (gfortran ≥ 10 иначе отвергает старые вызовы в стиле LAPACK);
4. сборка штатными makefile GACODE, **последовательно**: в них нет зависимостей по модулям Fortran, и `make -j` падает (`tglf_allocate.o` раньше `tglf_modules.o`). TGLF: `make -C tglf/src tglf_lib.a` (~10 с). NEO: `shared/{math,nclass,UMFPACK}`, `f2py/{geo,expro}` (именно в этом порядке: `expro` использует `geo.mod`), затем `neo/src neo_lib.a` (~5 с);
5. раскладка — та, что ждёт `exe/astra_rc` (каталог версии берётся из `astra_rc`, сейчас `jun25`/`dec24`):

   ```
   $ASTRA_EXT/tglf/jun25/lib/libtglf.a      (= tglf/src/tglf_lib.a)
   $ASTRA_EXT/tglf/jun25/inc/tglf*.mod
   $ASTRA_EXT/neo/dec24/lib/{neo,nclass,UMFPACK,math,expro,geo}_lib.a
   $ASTRA_EXT/neo/dec24/inc/*.mod
   ```

   Рядом файл `GACODE_COMMIT`: при совпадении коммита повторная сборка пропускается;
6. включение (выполняется при каждом запуске установщика, если библиотеки уже лежат в `$ASTRA_EXT`, даже без ключей):
   - в `exe/Makexpr` в цель `all:` возвращаются `$(XPR)/tglfi` и/или `$(XPR)/neo` (шаг 5 их убрал);
   - в `platform/env.ubuntu` дописывается `TGLF_FLAG=-DTGLF_INSTALLED`. Это обход дефекта upstream: `exe/Makesub` компилирует `sbr/tglf_serial.F90` с `$(TGLF_FLAG)`, но нигде его не определяет (в отличие от `RABBIT_FLAG`, `QLKNN_FLAG` и др.), поэтому `tglf_serial` всегда собирался заглушкой «not installed»;
   - там же `ASTRA_MAX_NWORKERS=${ASTRA_MAX_NWORKERS:-64}` и `MAX_NWORKERS=$ASTRA_MAX_NWORKERS`, а в `exe/Build` перед `echo 'Max nworkers'` вставляется строка `export MAX_NWORKERS=${ASTRA_MAX_NWORKERS:-$MAX_NWORKERS}`. Причина: обёртки `tglf_ipc`/`neo_ipc` жёстко запускают `nworkers=64` процессов и останавливают расчёт (`>>> ERROR nworkers= 64 exceeds #physical CPUs= 12`, `Quitting ASTRA`), если `MAX_NWORKERS < 64`. Без переменной обёртка берёт половину логических ядер, а `Build` в режиме с окном сам ставит число физических ядер. На ПК 64 процесса просто делят ядра.
7. проверка: после `flux_feqis` прогоняется `flux_tglf`/`aug34954`, t = 4–4.1 с (3 шага, 2 вызова TGLF; `TORBA`/`RABBIT` в этой модели уже закомментированы). Успех — `ASTRA normal exit`, строка `XPR wall time` (печатается после каждого вызова TGLF) и CDF. Только CDF недостаточно: при остановке по `MAX_NWORKERS` ASTRA успевает записать CDF с одним шагом.

Модели `flux_tglf` и `flux_neo` используют `IPEQL=5` (FEQIS) и закомментированные `TORBA`/`RABBIT`, поэтому идут без закрытых кодов. `flux_neo_tglf` требует SPIDER (`IPEQL=4`), `tglf_rot_pred_08` и другие — проверять по месту (`TORBA`/`RABBIT` закомментировать).

**NEO — экспериментально.** Библиотеки собираются, `xpr/neo` линкуется (проверено в 22.04 и 24.04), но расчёт `flux_neo` на стенде (24 логических ядра, 30 ГБ ОЗУ) за 11 минут не прошёл первый вызов: рабочие процессы `xpr/neo` занимали по 1.3–4.5 ГБ, система ушла в своп, расчёт остановлен вручную. NEO рассчитан на кластер (`exe/Build` при `sbatch` ставит `MAX_NWORKERS=128`).

## 4. Что делает `install_astra8.sh`

| шаг | действие | повторный запуск |
|---|---|---|
| 0 | проверка ОС, подтверждение | — |
| 1 | `apt-get install` пакетов из таблицы раздела 2 | apt пропускает установленное |
| 2 | по HTTPS с LF скачивается **только нужный коммит** (`git fetch --depth 1`, ~90 МБ вместо ~1 ГБ истории; `--full-history` — вся история) | `fetch` + сброс отслеживаемых файлов, затем патчи заново |
| 3 | сборка json-fortran в `$PREFIX/ext/json/9.0.2`; при `--blas mkl` — MKL с PyPI в `$PREFIX/ext/mkl/` | пропускается, если библиотека есть |
| 4 | `env.ubuntu`, `env.mypc`, правка `get_platform`, патч парсера | идемпотентно |
| 5 | копирование `repo_user_area/*` в рабочую копию, отключение xpr-целей и TORBA/RABBIT | `exe/` обновляется, в `equ/`, `exp/`, `udb/` ... дописываются только **отсутствующие** файлы |
| 5c | только с `--with-tglf`/`--with-neo`: OpenMPI из apt, GACODE (sparse, один коммит), сборка TGLF/NEO в `$PREFIX/ext/{tglf,neo}/`; если библиотеки есть — возврат целей в `Makexpr`, `TGLF_FLAG` и `MAX_NWORKERS` в `env.ubuntu`, строка в `exe/Build` (раздел 3.10) | сборка пропускается при том же коммите; включение повторяется и без ключей |
| 6 | команда `~/.local/bin/astra8`; если `~/.local/bin` нет в `PATH`, строка добавляется в `~/.bashrc` (стандартный `~/.profile` Ubuntu добавляет его только при следующем входе в систему) | перезаписывается, строка в `.bashrc` не дублируется |
| 7 | `make clean`, полная сборка и batch-расчёт `flux_feqis`/`aug34954`, t = 4–5 с; с `--with-tglf` ещё `flux_tglf`/`aug34954`, t = 4–4.1 с | выполняется заново (`--no-test` — пропустить) |

Параметры:

```
--prefix DIR        корень установки (по умолчанию ~/astra)
--ref REF           ветка, тег или ПОЛНЫЙ хэш коммита ASTRA (по умолчанию проверенный e66f6ba1...;
                    main — свежий код); короткий хэш работает только с --full-history
--full-history      скачать всю историю git ASTRA (~1 ГБ)
--blas openblas|mkl
--no-apt            не ставить пакеты
--no-test           не запускать тестовый расчёт
--with-tglf         собрать TGLF (GACODE), включить xpr/tglfi и tglf_serial
--with-neo          собрать NEO (GACODE), включить xpr/neo (экспериментально)
--gacode-ref REF    ветка, тег или ПОЛНЫЙ хэш GACODE (по умолчанию проверенный f34a9dd2...)
-y                  без вопросов
```

Журналы: `$PREFIX/install.log` (весь вывод apt/cmake/git, в том числе сборки GACODE), `$PREFIX/test_flux_feqis.log` (сборка и расчёт), с `--with-tglf` — `$PREFIX/test_flux_tglf.log`.

Изменения в upstream-файлах ASTRA (видны в `git -C ~/astra/a8 status`/`diff`):

- `get_platform` — 3 строки с `ASTRA_PLATFORM`;
- `pyparse/exp_parser.py` — патч из раздела 3.8;
- новые файлы `platform/env.ubuntu`, `platform/env.mypc`.

Всё остальное — в игнорируемой пользовательской области.

## 5. Проверка

### 5.1. Тестовый расчёт

`flux_feqis`/`aug34954`, t = 4–5 с, batch, без RABBIT/TORBEAM:

```
>>> ASTRA normal exit >>>
    Total time steps        21
    Average time step   50.000 msec
    >>> Astra wall time             00:00:30
```

Время на стенде (24 ядра, WSL2):

| этап | время |
|---|---|
| полная установка на чистой Ubuntu 24.04, включая apt (контейнер) | 2 мин 29 с |
| то же, Ubuntu 22.04 | 2 мин 52 с |
| скачивание ASTRA одним коммитом / с полной историей | 16 с / ~2.5 мин |
| `--blas mkl`: + скачивание MKL с PyPI | + ~1.5 мин |
| повторный запуск установщика с пересборкой и тестом | 42 с |
| повторный запуск `--no-test` | < 1 с |
| сборка ASTRA (`-O0`, 205 файлов) | ~10 с |
| тестовый расчёт `flux_feqis`, 21 шаг | 30 с |

### 5.2. Окно с графиками

`astra8 -m flux_feqis -v aug34954 -s 4 -e 5` открывает X11-окно «Per aspera ad ASTRA» (заголовок `ASTRA 8.6 -- Model: flux_feqis -- Data: aug34954`) с 8 профилями и стандартным меню (Graphic mode / Select / Control / In/Out / Status) — так же, как ASTRA 7.02 на машине лаборатории. В WSL2 на Windows 11 окно выводит WSLg, настраивать ничего не нужно.

### 5.3. Сравнение с эталоном

В `regressions/*.CDF` лежат эталонные результаты разработчиков (Intel + MKL на кластере IPP). Сравнение:

```bash
cd ~/astra/a8
python3 compareRegressions.py -m <model> -v <exp>   # отчёт: regressions/<exp><model>.log
```

Встроенный допуск `1e-7` рассчитан на ту же сборку, поэтому для gfortran он почти всегда показывает «ошибки». Смотреть надо на величину расхождения:

- **`fbe` / `AUG33040_2500`** (свободная граница, без закрытых модулей; `-s 2.48 -e 2.7`): `IPL`, `BTOR`, `TE`, `TI`, `NE` совпадают точно (заданы из эксперимента); `CU`, `MU` расходятся на ~0.1 %, `psiaxis` — на 1–3 % на всех шагах до последнего. На последнем шаге (t = 2.62) расчёт останавливается с `Error in find new boundary` — **так же, как в эталоне** (эталон тоже содержит 15 моментов времени, 2.48–2.62).
- **`flux_feqis` / `aug34954`**: с эталоном напрямую не сравнивается: эталон считан **с** RABBIT и TORBEAM (нагрев NBI и ECRH), а у нас они выключены.

- **OpenBLAS против MKL** (одна и та же сборка gfortran): на `fbe` 1013 из 1137 переменных совпадают побитово, остальные — до 1.5e-6; на `flux_feqis` — до 3e-6. Значит, расхождение с эталоном на уровне 1e-3 даёт компилятор (gfortran против Intel, разная арифметика с плавающей точкой), а не библиотека линейной алгебры.

### 5.4. Чистые системы

Установщик запускался без изменений в чистых контейнерах `ubuntu:24.04` и `ubuntu:22.04` (от root, `-y`): в обоих случаях код возврата 0, `ASTRA normal exit`, расчёт 31 с и 42 с. Отдельно проверены: обновление существующей полной копии (`~/astra` на стенде), `--ref main` и повторный запуск, сохраняющий файлы пользователя в `equ/`, `exp/`.

### 5.5. TGLF (`--with-tglf`), 2026-10-02

Чистые контейнеры, от root, `install_astra8.sh -y --with-tglf` (в 22.04 ещё `--with-neo`), ASTRA `e66f6ba1`, GACODE `f34a9dd2` (master 2026-09-30), стенд — 24 логических / 12 физических ядер, WSL2:

| | Ubuntu 24.04 | Ubuntu 22.04 |
|---|---|---|
| gfortran / OpenMPI | 13.3.0 / 4.1.6 | 11.4.0 / 4.1.2 |
| вся установка, включая apt и оба теста | 4 мин 0 с | 5 мин 28 с (с NEO) |
| сборка TGLF / NEO | ~10 с / — | ~10 с / ~5 с |
| `flux_feqis`, t = 4–5 с | 32 с | — |
| `flux_tglf`, t = 4–4.1 с, 3 шага | 1 мин 10 с | 2 мин 8 с |
| один вызов TGLF (64 процесса `xpr/tglfi`) | 28.7 с / 29.3 с | 55.2 с / 56.0 с |
| повторный запуск без ключей (`--no-apt --no-test`) | — | 20 с, TGLF и NEO остались включены, строки в `env.ubuntu`/`Build` не задвоились |

Результат `ncdf_out/aug34954flux_tglf.CDF`: `CAR17` = `tglf_out%chi_i`, `CAR18` = `chi_e` — конечные значения порядка 1–10 м²/с по радиусу (на последнем шаге, каждая 10-я точка: `6e-4 2.6 7.5 2.5 4.4 5.3 2.8 0.8 9.9 0`; ноль — на границе). Физически результаты не проверялись (эталона с TGLF для gfortran нет).

Последовательный вариант `flux_tglf_serial`, t = 4–4.1 с (24.04, ручная проверка той же раскладки): `ASTRA normal exit`, CDF записан, один вызов `tglf_serial` — 273 с (64 точки по ~4.3 с на одном ядре), всего 9 мин 14 с. Без `TGLF_FLAG` эта модель останавливалась бы на заглушке.

Режим с окном (GUI) с TGLF не запускался: правка `exe/Build` для него проверена только по тексту.

## 6. Docker

Папка `docker/` — альтернатива установке в систему. Образ собирается **тем же** `install_astra8.sh`, поэтому всё описанное выше относится и к нему.

```
docker/astra-docker.sh <аргументы astra8>      (на компьютере)
 ├─ docker build, если образа astra8 ещё нет (UID/GID = текущий пользователь)
 └─ docker run --rm
      -v ~/astra-data:/data                      модели, данные, результаты
      -e DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix  окно (WSLg / X11 рабочего стола)
      astra8 <аргументы>
         └─ entrypoint.sh  (в контейнере)
              ├─ /data/{equ,exp,udb,ncdf_out}: дописать примеры из образа, не перезаписывая
              ├─ ~/astra/a8/{equ,exp,udb,ncdf_out} -> симлинки на /data
              └─ astra8 <аргументы>   (или bash для "shell")
```

Устройство образа (`docker/Dockerfile`):

- база `ubuntu:${UBUNTU_VERSION}` (24.04; 22.04 проверен);
- пользователь `astra` с `sudo` без пароля. UID/GID задаются при сборке, чтобы файлы в `~/astra-data` принадлежали пользователю хоста. В `ubuntu:24.04` UID 1000 занят пользователем `ubuntu`, его удаляют;
- `install_astra8.sh -y --ref $ASTRA_REF --blas $BLAS` от имени `astra`. Внутри сборки проходит тестовый расчёт, так что неработающий образ просто не соберётся;
- после установки удаляются исходники json-fortran и списки apt.

Размер: ~2.2 ГБ по `docker images` (содержимое ~1.4 ГБ: `/usr` ~1.1 ГБ компиляторов и библиотек, ASTRA ~330 МБ). С полной историей git было бы ~4 ГБ. Сборка ~6 минут, в основном apt.

Проверено на стенде: сборка для 24.04 и 22.04; расчёт без окна с результатом в `~/astra-data/ncdf_out` (файлы принадлежат пользователю хоста); окно ASTRA из контейнера через WSLg; проброс произвольной команды (`astra-docker.sh ls /data`).

На Linux с рабочим столом (не WSL) `astra-docker.sh` вызывает `xhost +local:`, чтобы контейнер мог открыть окно на X-сервере. Это разрешает подключаться к X-серверу любым локальным процессам и контейнерам до перезапуска сеанса.

## 7. Ограничения и что дальше

- **Закрытые модули IPP** (RABBIT, TORBEAM, SPIDER, STRAHL): без них нет NBI-нагрева, ECRH и связки со SPIDER. Доступ — только через сотрудничество с IPP (аккаунт `gitlab.mpcdf.mpg.de`). Когда доступ есть, их можно собрать штатным `modules_install.sh` в тот же `$ASTRA_EXT` (раскладка каталогов — в `exe/astra_rc`).
- **TGLF** ставится ключом `--with-tglf` (раздел 3.10). **NEO** (`--with-neo`) собирается, но расчёт на ПК упирается в память (64 процесса по 1–4 ГБ). **QuaLiKiZ / QLKNN**: открытые, не добавлены (QuaLiKiZ с gitlab.com, нужен MPI, как NEO).
- **`showdata` + `readme`** (hello world из ASTRA 7): `exp/readme` в репозитории — файл формата A7 со ссылками на U-файлы `udb/00000.tes`, `udb/00000.tis`, которых в репозитории нет. После патча парсера он доходит до разбора ASCII-блоков профилей и падает на проверке размеров (`AssertionError` в `parse_exp2d`): формат блоков A7 отличается от A8. Нужны рабочие файлы лаборатории из установки ASTRA 7.
- **Оптимизация**: `FC_FLAGS` содержит `-O0`, как у разработчиков (они используют его и для Intel). Ускорение через `-O2` возможно, но не проверялось на совпадение результатов.

## 8. Частые ошибки

| симптом | причина | что делать |
|---|---|---|
| `cannot execute: required file not found` | CRLF в скриптах | переклонировать с `core.autocrlf=false` (установщик делает сам) |
| `platform/env.xxx: No such file or directory` | hostname совпал с чужой площадкой или нет env-файла | запускать через `astra8` (ставит `ASTRA_PLATFORM=ubuntu`) |
| `Invalid form of array reference` | нет `-ffree-line-length-none` | проверить `FC_FLAGS` в `platform/env.ubuntu` |
| `undefined reference to XOpenDisplay` | `--as-needed` | проверить `LD_FLAGS` в `platform/env.ubuntu` |
| `libmkl_core.a: No such file` | выбран MKL, но он не скачался | переустановить с `--blas openblas` или повторить `--blas mkl` |
| `Error from exe/Makexpr-compilation` | цели TGLF/QLK/NEO без библиотек | перезапустить установщик (уберёт цели) |
| `>>> ERROR nworkers= 64 exceeds #physical CPUs=...`, `Quitting ASTRA` | модель с `TGLF_IPC`/`NEO_IPC`, а `MAX_NWORKERS` < 64 | поставить `--with-tglf` (пишет `ASTRA_MAX_NWORKERS=64` в `env.ubuntu`) или `export ASTRA_MAX_NWORKERS=64` |
| `TGLF not installed`, `TGLF_SERIAL` останавливается | нет `$ASTRA_EXT/tglf/jun25/lib/libtglf.a` или `TGLF_FLAG` | `bash install_astra8.sh --with-tglf` |
| `cannot fetch GACODE '...'` | короткий хэш в `--gacode-ref` или нет доступа к GitHub | полный хэш (40 символов); проверить `git ls-remote https://github.com/gafusion/gacode.git` |
| `FileNotFoundError: .../udb/...` | файл эксперимента ссылается на отсутствующий U-файл | положить U-файл в `~/astra/a8/udb/` |
| окно не появляется, `cannot open display` | нет X-сервера | Windows 11 + WSL2 (WSLg) или рабочий стол Ubuntu; в Windows 10 нужен X-сервер (VcXsrv) и `export DISPLAY=...`; или считать с `-batch` |
| расчёт «молча» не идёт, код возврата 0 | `as_exe` не передаёт ошибки | смотреть вывод: ищите `Error` выше по тексту |
| `cannot fetch '...'` при установке | короткий хэш коммита в `--ref` | дать полный хэш (40 символов) или добавить `--full-history` |
| `astra8: command not found` сразу после установки | терминал открыт до установки | открыть новый терминал |
| Docker: `docker: command not found` в WSL | не включена интеграция Docker Desktop с этой Ubuntu | Docker Desktop → Settings → Resources → WSL integration |
| Docker: `Permission denied` в `~/astra-data` | образ собран с другим UID | пересобрать через `astra-docker.sh` (подставляет ваш UID) или с `--build-arg USER_UID=$(id -u)` |
