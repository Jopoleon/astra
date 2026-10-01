# Площадки ASTRA 8: `get_platform` и `platform/env.*`

Справочник по тому, как ASTRA 8 подбирает настройки компилятора и библиотек под конкретный компьютер, и что лежит в 22 файлах `platform/env.*` upstream-репозитория.

Источник: `astra-src/`, ветка `local` (по этим файлам совпадает с `main`, коммит `6728111`, 2026-10-01). Все пути ниже указаны относительно `astra-src/`, если не оговорено иное. Где вывод сделан по косвенным признакам, стоит пометка «предположительно».

## 1. Как ASTRA выбирает площадку

### 1.1. `get_platform`

`get_platform` — bash-скрипт в корне рабочей копии, который печатает имя площадки. Используются три признака:

- `AWD` — абсолютный путь к каталогу, где лежит сам скрипт (`dirname $0`);
- `node` — `uname -n` (hostname);
- `os` — `uname`.

Правила проверяются по порядку, срабатывает первое совпадение:

| # | условие | результат | есть `env.*`? |
|---|---|---|---|
| 1 | `AWD` начинается с `/localhome` | `goblin` | **нет** |
| 2 | `AWD` начинается с `/home/ITER` | `iter` | да |
| 3 | `AWD` начинается с `/home/matlab` | `docker` | да |
| 4 | `AWD` начинается с `/global/homes` | `perlmutter` | да |
| 5 | `AWD` начинается с `/compass` | `cz` | да |
| 6 | hostname `*mit.edu` **или** `node*` | `mit` | да |
| 7 | hostname `omega*` | `omega` | да |
| 8 | hostname `cx-ld-prod*` или `hz-ld-*` | `aug` | да |
| 9 | hostname `viz*` | `gway` | да |
| 10 | hostname `tok*` | `tok` | да |
| 11 | hostname `iris*` | `iris` | **нет** (удалён 2025-08-07, `7def936b`) |
| 12 | hostname `ncku-first` | `ncku` | да |
| 13 | hostname `smartssv*` | `sevilla` | да |
| 14 | hostname `rat2*` | `rat2` | да |
| 15 | hostname `puhti*` | `puhti` | да |
| 16 | hostname `freia*` | `freia` | да |
| 17 | hostname `fc-deb*` | `hgw` | да |
| 18 | hostname `feynman` | `columbia` | да |
| 19 | hostname `tohtorisms` | `tohtori` | да |
| 20 | `uname` = `Darwin` | `darwin` | да |
| — | иначе | `mypc` | **нет** |

Затем, если в окружении **задана** (пусть даже пустая) переменная `ASTRA_COMPILER`, к имени добавляется суффикс: `platform=${platform}_${ASTRA_COMPILER}` (проверка `[[ -v ASTRA_COMPILER ]]`). Переменную выставляет `exe/as_exe.py` из ключа `-c/--compiler` (`gcc`, `ifx`); значение запоминается в `tmp/astra_log` и при следующем запуске без `-c` берётся оттуда (`repo_user_area/exe/as_exe.py`, строки 68, 85–115). Варианты с суффиксом есть только для TOK: `env.tok_gcc`, `env.tok_ifx`.

Файлы `env.lac` и `env.w7x` из `get_platform` **недостижимы**: правила для них нет (правило `fc-deb*` до 2025-03-05 вело на `w7x`, теперь на `hgw`, коммит `4cd56f8a`). Использовать их можно, только поправив `get_platform` или выставив переменную-переключатель, как это делает наш установщик (см. 3.3).

### 1.2. Если ни одно правило не сработало

Скрипт возвращает `mypc`, а `platform/env.mypc` в репозитории нет (и не было в истории). То же для `goblin` и `iris`. Дальше:

- `source $AWD/platform/env.mypc` в bash не фатален: печатается `No such file or directory`, скрипт продолжает работу;
- `FC`, `FC_FLAGS`, `MKL_LIBDIR`, `ASTRA_EXT` и т.д. остаются пустыми; `make` берёт свой встроенный `FC` (`f77`), сборка падает на компиляции/линковке;
- отдельный риск — ложные совпадения: любой хост с именем `node*` считается MIT, `tok*` — IPP TOK, `omega*` — GA, `viz*` — Gateway. Подтягиваются чужие абсолютные пути. Подробнее: `installer/docs/DETAILED.md`, 3.2.

`modules_inst_loc.sh` (установка внешних модулей на ПК, коммит `7698c2f7`, 2026-08-10) содержит ветку `if [[ "$platform" == "linux" ]]` с `apt-get`, но `get_platform` никогда не выдаёт `linux`, а `env.linux` нет. Предположительно, автор у себя добавляет такой файл и правит `get_platform` локально; в upstream Linux-ветка не срабатывает.

### 1.3. Кто читает `env.<platform>`

Файл подключается через `source`, то есть выполняется в текущей оболочке, со всеми побочными эффектами (`module load`, `mkdir`, `echo`...).

| скрипт | что делает после `source` |
|---|---|
| `install.sh` | берёт `ASTRA_EXT`, печатает его, `make -f exe/Makefile clean`, тестовый прогон `flux_feqis`/`aug34954` |
| `repo_user_area/exe/as_exe` | добавляет `PYTHON_BINDIR` в `PATH`, запускает `as_exe.py` |
| `repo_user_area/exe/Build` | затем `source exe/astra_rc` (или `exe/astra_rc_${ASTRA_COMPILER}`), `ulimit -s`, `make -f exe/Makefile`, `make -f exe/Makexpr`, запуск или `sbatch` |
| `compReg.sh` | прогон регрессий; проверяет `PYTHON_BIN`, но в `PATH` добавляет `PYTHON_BINDIR` (опечатка) |
| `modules_install.sh` | сборка внешних модулей в `$ASTRA_EXT` на кластерах; принимает аргумент `comp` и **сам** добавляет суффикс `_${comp}` |
| `modules_inst_loc.sh` | то же для ПК (apt/brew); суффикс добавляет, только если файл `env.${platform}_${comp}` существует |

Кроме того, `exe/Makefile` и `exe/Makexpr` сами вызывают `get_platform` (`PLATFORM := $(shell $(AWD)/get_platform)`), чтобы выбрать OpenBLAS вместо MKL и путь к X11 при `PLATFORM=darwin`.

При запуске одна и та же модель проходит `env.*` дважды: в `as_exe` и в `Build`.

Каталог `exe/` в рабочей копии — это копия `repo_user_area/exe/`, которую делает `install.sh`.

### 1.4. Контракт переменных

Что должен (или может) определить `env.<platform>` и кто это потребляет.

| переменная | обязат. | кто использует | назначение |
|---|---|---|---|
| `FC` | да | `exe/Makefile`, `Makeastra`, `Makesub`, `Makexpr` | компилятор Fortran. В `Makexpr` сравнивается **строго** с `gfortran`: иначе считается Intel (`-qopenmp`, `-lifcore`, линковка через `$(CC)`) |
| `CC` | да | `Makeastra`, `Makesub`, `Makexpr` | компилятор C (`libgraph`, IPC-обёртки) |
| `MPIFC` | для `xpr/` | `Makexpr` (NEO, QuaLiKiZ), `modules_install.sh` (QuaLiKiZ, GACODE) | MPI-обёртка Fortran |
| `MPICC` | для `xpr/` (Intel) | `Makexpr`: `LD_MPI := $(MPICC) -qopenmp` | линковщик IPC-программ в Intel-ветке |
| `FC_FLAGS` | да | все `Make*` | флаги компиляции; `Makefile` дописывает `-I$(AWD)` |
| `LD_FLAGS` | да | `exe/Makefile`, линковка `bin/*.exe` | флаги линковки (`-I$(AWD)` дописывается) |
| `MKL_LIBDIR` | да, кроме `darwin` | `Makefile`, `Makexpr` | каталог со **статическими** `libmkl_intel_lp64.a`, `libmkl_core.a`, `libmkl_sequential.a` |
| `OPENBLAS_LIBDIR` | только `darwin` | `Makefile`, `Makexpr` | `-L... -lopenblas` вместо MKL |
| `COMP_LIBDIR` | для Intel | `Makexpr` (`-lifcoremt -lifport -lifcore`), `astra_rc` (`IRC_LIB=-lirc` для `libtorbeamB.so`) | runtime-библиотеки Intel |
| `ASTRA_EXT` | да | `astra_rc*`, `install.sh`, `modules_*install*.sh`; `src/astra/read_input.f90:68` (`getenv`, далее, предположительно, не используется) | корень внешних модулей |
| `LD_LIBRARY_PATH` | обычно | runtime; `Build` переносит его в `batch.sh` | почти все Intel-конфиги **перезаписывают** его значением `COMP_LIBDIR` |
| `PATH` | обычно | — | каталоги `bin` компилятора и MPI (`COMP_BINDIR`, `MPI_BINDIR`) |
| `PYTHON_BINDIR` | желат. | `exe/as_exe`, `compReg.sh` | каталог с `python3` (numpy, scipy, netCDF4) для парсера и `json2cdf` |
| `PYTHON_BIN` | опц. | `src/astra/a2vmec.f90:41` | интерпретатор для связки с VMEC/STELLOPT (стеллараторная ветка) |
| `STELLOPT_HOME` | опц. | `a2vmec.f90:42`, `vmec/python/*`; по умолчанию `astra_rc`: `${ASTRA_EXT}/STELLOPT/` | STELLOPT |
| `MPI_COMMAND` | опц. | `a2vmec.f90:40` | команда `mpirun` для VMEC |
| `SL_QUEUE` | опц. | `exe/Build` | **дословно** вставляется строкой в `exe/batch.sh` при `-batch` и наличии `sbatch`; с 2026-02-04 (`055b5a46`) это должен быть готовый блок `#SBATCH ...` |
| `I_MPI_ROOT`, `I_MPI_FABRICS` | нет | только Intel MPI runtime | `I_MPI_FABRICS=shm` — только разделяемая память, т.е. один узел |
| `MKLROOT` | нет | `modules_install.sh` (GACODE, gfortran) | |
| `NETCDF_HOME`, `NETCDF_LIB` | нет | только закомментированная строка `exe/Makefile:157` | в основной сборке ASTRA **netCDF не линкуется**; нужен RABBIT и STRAHL, их собирают `modules_*install*.sh` |

Переменные, которые задаёт **не** `env.*`, а `exe/astra_rc` по содержимому `$ASTRA_EXT` (см. 1.5): `JSON_LIB`, `JSON_INC`, `SPIDER_LIB/INC`, `TORBEAM_LIB`, `RABBIT_LIB/INC`, `TORIC_LIB/INC`, `QLK_LIB/INC`, `QLKNN_LIB/INC`, `TGLF_LIB/INC`, `NEO_LIB/INC`. `Build` экспортирует `N_NODES`, `AWD`, `MAX_NWORKERS`, `X11`.

Объявлены в некоторых `env.*`, но ничем не читаются: `NODES`, `CPUS`, `TIME` (`env.tohtori`), `PYTHONHOME` (читает сам Python, `env.gway`), `OPENMPI_LIB`, `LIBOMP_LIB`, `OPENBLAS_INCLUDEDIR` (`env.darwin`), `INTEL_LICENSE_FILE` (читает компилятор, `env.cz`).

### 1.5. `exe/astra_rc`: где ищутся внешние коды

`repo_user_area/exe/astra_rc` одинаков для всех площадок: он не определяет, что «есть» на площадке, а проверяет наличие файлов в `$ASTRA_EXT` и при отсутствии делает `unset` (с коммита `2ec99235`, 2026-05-12, ASTRA собирается с заглушками).

| код | путь (по умолчанию, Intel) | `astra_rc_gcc` | проверка наличия |
|---|---|---|---|
| json-fortran | `json/9.0.2` | `json/gcc_9.0.2` | **нет**, библиотека обязательна |
| SPIDER | `spider/latest` | `spider/gcc_may26` | да |
| TORBEAM | `torbeam/nov24/lib/libtorbeamB.so` | `torbeam/gcc_apr26` | да |
| RABBIT | `rabbit/may26/lib/librabbit.{so,dylib}` | `rabbit/gcc_may26` | да |
| TORIC | `toric/apr26/intel/2025.1` + FFTW из `/mpcdf/soft/SLE_15/packages/znver4/fftw/...` | `toric/apr26/gfortran/2025.1` | да |
| QuaLiKiZ | `qualikiz/nov24`, `libQLK-intel-release-default.a` (или `-gcc-`) | `qualikiz/gcc_apr26` | **нет** |
| QLKNN | `qlk_nn/nov24` | `qlk_nn/gcc_apr26` | да |
| TGLF | `tglf/jun25/lib/libtglf.a` | `tglf/gcc_apr26` | да |
| NEO | `neo/dec24/lib/*.a` | `neo/gcc_apr26` | **нет** |
| STELLOPT | `${STELLOPT_HOME:-$ASTRA_EXT/STELLOPT/}` | не задаётся | — |

На `darwin` `astra_rc` берёт **безверсионные** каталоги (`$ASTRA_EXT/json`, `/rabbit`, ...) — именно так их раскладывает `modules_inst_loc.sh`. `astra_rc_ifx` фиксирует те же версии, что и основной, без автоопределения gcc-библиотек.

Поэтому колонка «внешние коды» в таблице раздела 2 говорит только о том, что следует из самого `env.*`. Реальный набор определяется содержимым `$ASTRA_EXT` на машине, а его из репозитория не видно. `modules_install.sh` копирует использованный `env.<platform>` внутрь каждого собранного модуля (`cp $AWD/platform/env.$platform $RABBIT_INSTALL/` и т.п.), поэтому на площадке можно проверить, каким конфигом модуль собран.

## 2. Сводная таблица `platform/env.*`

### 2.1. Где это

| файл | как распознаётся | организация / машина | последняя правка |
|---|---|---|---|
| `env.aug` | hostname `cx-ld-prod*`, `hz-ld-*` | IPP Garching, Linux-кластер группы ASDEX Upgrade (`/shares/departments/AUG`) | 2026-05-13 |
| `env.tok` | hostname `tok*` | IPP / MPCDF, кластер TOK (SLES 15, `znver4`, Slurm `p.tok*`). Основная площадка разработчиков | 2026-09-18 |
| `env.tok_ifx` | `tok*` + `ASTRA_COMPILER=ifx` | то же, явный вариант Intel ifx | 2026-05-13 |
| `env.tok_gcc` | `tok*` + `ASTRA_COMPILER=gcc` | то же, GNU | 2026-05-13 |
| `env.hgw` | hostname `fc-deb*` | IPP Greifswald (W7-X), пользователь `papo` | 2026-05-13 |
| `env.w7x` | **недостижим** | IPP Greifswald, личная установка AStar3D (`$HOME/ipp-hgw/AStar3D`), предположительно ветка `stars3d` | 2026-08-05 |
| `env.iter` | путь `/home/ITER*` | ITER Organization, кластер SDCC (EasyBuild `/work/imas`) | 2026-05-13 |
| `env.gway` | hostname `viz*` | EUROfusion Gateway (`/afs/eufus.eu`, модули CINECA) | 2026-05-13 |
| `env.perlmutter` | путь `/global/homes*` | NERSC Perlmutter, проект `m2116` | 2026-05-13 |
| `env.omega` | hostname `omega*` | General Atomics, кластер omega (`/fusion/usc`, `/fusion/projects`) | 2026-05-13 |
| `env.mit` | hostname `*mit.edu` или `node*` | MIT PSFC, кластер ORCD/Engaging (`/orcd`, очередь `sched_mit_psfc_r8`) | 2026-05-13 |
| `env.columbia` | hostname `feynman` | предположительно Columbia University (машина `feynman`, `/mnt/codes`) | 2026-05-19 |
| `env.cz` | путь `/compass*` | IPP Prague (Czech Academy of Sciences), токамак COMPASS | 2026-05-13 |
| `env.rat2` | hostname `rat2*` | предположительно Consorzio RFX, Падуя (кластер `rat2`, RHEL 7) | 2026-05-13 |
| `env.puhti` | hostname `puhti*` | CSC (Финляндия), суперкомпьютер Puhti, проект `project_2007159` (предположительно группа VTT) | 2026-04-16 |
| `env.tohtori` | hostname `tohtorisms` | Финляндия; в коммите `f7c8cc7f` от 2025-10-22: «VTU Finnland», предположительно VTT | 2026-05-13 |
| `env.freia` | hostname `freia*` | UKAEA, Culham, кластер Freia | 2026-04-16 |
| `env.ncku` | hostname `ncku-first` | National Cheng Kung University (Тайвань), пользователь `astrainstaller` | 2026-05-13 |
| `env.sevilla` | hostname `smartssv*` | University of Seville, токамак SMART | 2026-05-13 |
| `env.lac` | **недостижим** | Swiss Plasma Center (EPFL), кластер `lac` (lac5, lac10), токамак TCV; пользователь `cooseman` | 2026-07-13 |
| `env.docker` | путь `/home/matlab*` | Docker-образ; предположительно образ с MATLAB-пользователем, ASTRA_EXT указывает на AUG | 2026-05-13 |
| `env.darwin` | `uname` = `Darwin` | macOS + Homebrew | 2026-07-30 |
| (нет файла) `goblin` | путь `/localhome*` | неизвестно, правило добавлено 2023-07-06 (`9b96fd57`) | — |
| (нет файла) `iris` | hostname `iris*` | GA, кластер iris; конфиг удалён 2025-08-07 | — |

Почти все правки — Giovanni Tardini (IPP); `env.tok`/`env.w7x` также правил Emiliano Fable (2026-08-01).

### 2.2. Чем собирают

Общие флаги Intel-конфигов: `FC_FLAGS="-r8 -no-prec-div -traceback -O0 -w -fPIC -qopenmp"`, `LD_FLAGS="-qopenmp -traceback"` — далее «стандарт Intel». У GNU: `FC_FLAGS="-O0 -w -fPIC -fopenmp"`, `LD_FLAGS="-fopenmp"` — «стандарт GNU».

| файл | FC / CC | MPIFC / MPICC | MPI | BLAS/LAPACK | netCDF | Python | `module load` | ASTRA_EXT | флаги |
|---|---|---|---|---|---|---|---|---|---|
| `aug` | ifort / icc (oneAPI 2023.1) | mpiifort / mpiicc | Intel MPI 2021.9 | MKL 2023.1 | — | miniforge3 2024.09 | нет, прямые пути `/shares/software/aug-dv/...` | `/shares/departments/AUG/users/git/ASTRA_LIBRARIES_EXT` | стандарт Intel |
| `tok` | ifx / icx (oneAPI 2025.1) | mpiifx / mpiicx | Intel MPI 2021.15; `MPI_COMMAND` — OpenMPI 5.0.8 | MKL 2025.1 (`MKLROOT`) | MPCDF `netcdf-mpi` 4.9.2 (`NETCDF_HOME`, `NETCDF_LIB`) | waterboa 2025.06 (`PYTHON_BINDIR`, `PYTHON_BIN`) | нет, пути `/tok/u/system/soft`, `/mpcdf/soft` | `/tokp/work/software/TOK_2023/software/ASTRA_LIBRARIES_EXT` | стандарт Intel; `SL_QUEUE`, `STELLOPT_HOME` |
| `tok_ifx` | ifx / icx (2025.1) | mpiifx / mpiicx | Intel MPI 2021.15 | MKL 2025.1 | — | waterboa 2024.06 | нет | как `tok` | стандарт Intel; `SL_QUEUE` |
| `tok_gcc` | gfortran / gcc 14 | mpifort / mpicc | Intel MPI 2021.15 (сборка под gcc 14) | MKL 2025.1 из модуля | — | waterboa 2024.06 | `module purge; module load gcc/14 mkl/2025.1` | как `tok` | стандарт GNU; `SL_QUEUE` |
| `hgw` | ifx / icx (2025.1) | mpiifx / mpiicx | Intel MPI 2021.15 | MKL 2025.1 | — | личный venv | нет, `/opt/intel/oneapi` | `~/ASTRA_LIBRARIES_EXT` | стандарт Intel |
| `w7x` | ifx / icx (2025.2) | **mpif90** / mpiicx | Intel MPI 2021.16; `MPI_COMMAND=/usr/bin/mpirun -n` | MKL 2025.2 | `/usr` + `$ASTRA_EXT/netcdf-parallel`, `hdf5` | `$HOME/myenv` (`PYTHON_BIN`) | нет, Intel лежит в `$HOME` | `$HOME/ipp-hgw/AStar3D` | **`FC_FLAGS`/`LD_FLAGS` не заданы** |
| `iter` | ifx / icx (2024.2) | mpiifx / mpiicx | Intel MPI 2021.13 | MKL 2024.2 (EasyBuild `imkl`) | — | `SciPy-bundle/2023.11` (Python 3.11.5) | `source /etc/profile.d/modules.sh; module load SciPy-bundle/...` | `/home/ITER/tarding/ASTRA_LIBRARIES_EXT` | стандарт Intel |
| `gway` | ifx / icx (2024.1) | mpiifx / mpiicx | Intel MPI 2021.12 | MKL 2024.0 | модули `netcdf-c/4.9.2`, `netcdf-fortran/4.6.1` | miniforge3 (`PYTHONHOME`) | `cineca-modules/1.0 intel-oneapi-compilers/2024.1.0 intel-oneapi-mkl/... netcdf-c/... netcdf-fortran/...` | `/afs/eufus.eu/user/g/g2gtardi/ASTRA_LIBRARIES_EXT` | стандарт Intel |
| `perlmutter` | ifort / icc (2022.1) | mpiifort / mpiicc | Intel MPI 2021.10 | MKL 2022.1 | — | NERSC conda 24.1.0, Python 3.11 | нет, `/opt/intel/oneapi` (не Cray PE) | `/global/cfs/cdirs/m2116/ASTRA_LIBRARIES_EXT` | стандарт Intel |
| `omega` | ifort / icc (Intel 2020) | mpifort / mpicc | Intel MPI 2021 U3 **из CST Studio Suite 2022** | MKL (Intel 2020) | — | GA conda `general` | нет | `/fusion/projects/codes/astra/c8/8.2/intel2020/ASTRA_LIBRARIES_EXT` | стандарт Intel |
| `mit` | ifx / icx (2024.2.1) | mpiifx / mpiicx | Intel MPI 2021.13 | MKL 2024.2 | — | anaconda3 2023.07 | нет, `/orcd/software/community/...` | `/orcd/pool/003/pablorf_shared/code_repo/ASTRA_LIBRARIES_EXT` | стандарт Intel; `SL_QUEUE` (старый формат) |
| `columbia` | ifx / icx (2025.3) | mpiifx / mpiicx | Intel MPI 2021.18 | MKL 2025.3 | — | `/usr/bin` | нет | `/mnt/codes/astra/ASTRA_LIBRARIES_EXT` | стандарт Intel |
| `cz` | ifort / icc (Intel 2018.3) | mpiifort / mpiicc | Intel MPI 2018 | MKL 2018 | — | miniforge3 (путь сломан, см. 4) | нет | `/compass/home/tardini/astra_ext` | стандарт Intel; `INTEL_LICENSE_FILE` |
| `rat2` | ifort / icc (Intel 2018u1) | mpiifort / mpiicc | Intel MPI 2018 | MKL 2018 | самосборный, `/home/tardini/soft/netCDF/build` в `LD_LIBRARY_PATH` | anaconda 2019.03, Python 3.7.1 | нет | `/home/tardini/ASTRA_LIBRARIES_EXT` | стандарт Intel |
| `puhti` | ifort / icc | mpifort / mpicc | OpenMPI (модуль) | MKL 2022.1 (жёсткий путь в spack) | — | модуль `python-data` | `module purge; module load python-data intel-oneapi-compilers intel-oneapi-mkl openmpi` | `/projappl/project_2007159/ASTRA_LIBRARIES_EXT` | стандарт Intel |
| `tohtori` | ifort / icc (бинарники 2019.4, библиотеки 2023.2) | mpiifort / mpiicc | Intel MPI 2021.10 | MKL 2023.2 | — | модуль `python/3.12` | `module purge; module load python/3.12` | `/home/akadam/builds/A88_EXT` | стандарт Intel; `SL_QUEUE`, `NODES`, `CPUS`, `TIME` |
| `freia` | ifort / icc (2020.3) | mpifort / **mpicc** (опечатка, см. 4) | OpenMPI 4.1.0 | MKL из модуля ifort (`$MKLROOT/lib/intel64_lin`) | — | `python/3.9-minibundle` | `module purge; module load ifort/2020.3.110 openmpi/4.1.0 openmpi-devel/4.1.0; module use ...; module load python/3.9-minibundle` | `$HOME/ASTRA_LIBRARIES_EXT` | стандарт Intel |
| `ncku` | ifort / icc (oneAPI 2024.2) | mpiifort / mpiicc | Intel MPI 2024.2 | MKL 2024.2 | — | не задан | нет | `/home/astrainstaller/ASTRA_LIBRARIES_EXT` | стандарт Intel |
| `sevilla` | `ifort -diag-disable=10448` / `icc ...` (2024.2) | `mpiifort ...` / `mpiicc ...` | Intel MPI 2024.2 | MKL 2024.2 | самосборный `/home/alfrodgon/netCDF` (как rpath, см. 4) | не задан | нет | `/home/software/ASTRA/ASTRA_EXT_NEW` | стандарт Intel; подавление предупреждения 10448 |
| `lac` | ifx / icx (oneAPI `latest`) | mpiifx / — | Intel MPI (неявно) | MKL `latest` | — | venv (`PYTHON_BIN` только) | нет | `/lacdata/users/cooseman/soft` | **`FC_FLAGS`/`LD_FLAGS` не заданы**; создаёт/удаляет каталоги |
| `docker` | ifort / `gcc -O` | mpiifort / mpiicc | Intel MPI `latest` | MKL `latest` | — | не задан | нет | `/shares/departments/AUG/...` (путь AUG) | стандарт Intel |
| `darwin` | `gfortran-N` / `gcc-N` из Homebrew (старшая версия) | mpif90 (или mpifort) / mpicc | Homebrew `open-mpi` | **OpenBLAS** (Homebrew) | Homebrew (ставит `modules_inst_loc.sh`) | `python` из `PATH` | нет (Homebrew, XQuartz) | `${AWD}/astra_modules` (внутри рабочей копии) | стандарт GNU |

Итого по компиляторам: ifort — 11 конфигов (aug, perlmutter, omega, cz, rat2, puhti, tohtori, freia, ncku, sevilla, docker), ifx — 9 (tok, tok_ifx, hgw, w7x, iter, gway, mit, columbia, lac), gfortran — 2 (tok_gcc, darwin). **nvfortran, flang, Cray `ftn` не используются нигде**, в том числе на Perlmutter.

### 2.3. Внешние коды, видимые из `env.*`

| что | где |
|---|---|
| RABBIT, TORBEAM, SPIDER, TGLF, NEO, QuaLiKiZ, QLKNN, json-fortran | ни в одном `env.*` не упоминаются; ищутся `exe/astra_rc` в `$ASTRA_EXT` (раздел 1.5) |
| TORIC | только TOK: путь FFTW в `astra_rc` жёстко указывает на MPCDF `znver4` |
| STELLOPT/VMEC (стеллараторная ветка `a2vmec`) | `tok` (`STELLOPT_HOME`, `NETCDF_*`, `MPI_COMMAND`, `PYTHON_BIN`), `w7x` (`NETCDF_*`, `MPI_COMMAND`, `PYTHON_BIN`, hdf5 и parallel netCDF в `$ASTRA_EXT`) |
| TCV-данные | `lac`: создаёт `$AWD/dat/tcvsimdat/` |
| RABBIT с `-xHost` | `tohtori`: комментарий о том, что из-за этого работают только узлы `gen04_epyc`, `gen05_skylake` |

Версии внешних модулей «по площадкам» косвенно видны в `ASTRA_EXT` omega: `.../astra/c8/8.2/intel2020/...` — предположительно сборка под ASTRA 8.2.

## 3. Общие закономерности и отличия

### 3.1. Типовой Intel-конфиг

20 Intel-конфигов из 22 — вариации одного шаблона (у `puhti`, `freia`, `gway`, `lac`, `w7x` часть блоков отсутствует) (сформирован коммитами «Consistent platform/env.*» 2025-08-07 и «Updated platform/env.*» 2026-04-16):

```bash
PKG=<корень Intel>
export I_MPI_ROOT=${PKG}/mpi/<ver>
export I_MPI_FABRICS=shm
export MKL_LIBDIR=${PKG}/mkl/<ver>/lib/intel64
export COMP_LIBDIR=${PKG}/compiler/<ver>/lib
export LD_LIBRARY_PATH=${COMP_LIBDIR}
export FC=ifx CC=icx MPIFC=mpiifx MPICC=mpiicx      # или ifort/icc/mpiifort/mpiicc
export FC_FLAGS="-r8 -no-prec-div -traceback -O0 -w -fPIC -qopenmp"
export LD_FLAGS="-qopenmp -traceback"
COMP_BINDIR=...; MPI_BINDIR=$I_MPI_ROOT/bin
[[ ":$PATH:" != *":$COMP_BINDIR:"* ]] && export PATH=${COMP_BINDIR}:${MPI_BINDIR}:${PATH}
export PYTHON_BINDIR=...
export ASTRA_EXT=...
```

Отличия между площадками — почти только пути и версии Intel. Систем модулей используют 6 конфигов (`tok_gcc`, `iter`, `gway`, `puhti`, `tohtori`, `freia`); остальные прописывают абсолютные пути вручную. Это осознанный выбор: коммит `f2c702b9` (2024-11-18) называется «Move away from modules».

### 3.2. Что одинаково везде

- **`-O0`** во всех конфигах: ASTRA и в продакшене собирается без оптимизации (оптимизируются только внешние модули, например GACODE с `-Ofast`).
- `-w` — все предупреждения подавлены.
- `-r8` (REAL по умолчанию 8 байт) — только в Intel-конфигах. GNU-конфиги (`tok_gcc`, `darwin`) **не** ставят `-fdefault-real-8`; по коммиту `c7202c76` (2026-04-20, «double precision everywhere») код переведён на явную двойную точность, так что `-r8` для ядра, видимо, уже не нужен. Наш `env.ubuntu` тоже без него, и регрессии совпадают (`installer/docs/DETAILED.md`, 5.3).
- MKL всегда **статический и последовательный** (`libmkl_sequential.a`) — параллелизм делает сама ASTRA через OpenMP.
- `I_MPI_FABRICS=shm` — Intel MPI ограничен одним узлом; IPC-программы `xpr/` (TGLF/NEO/QuaLiKiZ) работают внутри узла.
- netCDF основной сборке не нужен; в `env.*` он появляется только там, где самостоятельно собирают RABBIT или работает связка с VMEC.

### 3.3. Ближайший к обычной Ubuntu/gfortran конфиг

| кандидат | плюсы | минусы |
|---|---|---|
| `env.tok_gcc` | Linux, gfortran, `mpifort`/`mpicc`, стандарт GNU, сами разработчики проверяют его на TOK («gfortran working, with RABBIT and TORBEAM», 2026-04-16); с ним согласован `astra_rc_gcc` | MKL и Intel MPI из модулей MPCDF; `SL_QUEUE`; не хватает `-ffree-line-length-none` и `-Wl,--no-as-needed` |
| `env.darwin` | gfortran + OpenBLAS + OpenMPI из пакетного менеджера — по духу ровно то, что даёт apt | macOS-специфичен; OpenBLAS включается в `Makefile`/`Makexpr` только при `PLATFORM=darwin`; `FC=gfortran-N` ломает проверку `ifeq ($(FC),gfortran)` в `Makexpr` |

**Лучший шаблон — `env.tok_gcc`** (набор переменных и флаги), с подходом `env.darwin` к BLAS. Ровно так устроен наш `platform/env.ubuntu` (`installer/docs/DETAILED.md`, 3.2): переменные `tok_gcc` без `module`/Intel MPI, плюс `-ffree-line-length-none`, `-Wl,--no-as-needed`, `MAKEFLAGS="MKL_LIB=-lopenblas"` и `OPENBLAS_NUM_THREADS=1`. Выбор площадки установщик делает через переменную `ASTRA_PLATFORM`, которую проверяют три строки, вставленные в начало `get_platform`. При заданной `ASTRA_PLATFORM` скрипт выходит сразу, поэтому суффикс `ASTRA_COMPILER` к ней не добавляется.

## 4. Находки и странности

### 4.1. Мёртвые и недостижимые записи

| что | суть |
|---|---|
| `goblin`, `iris`, `mypc` | `get_platform` выдаёт имена, для которых нет `env.*`; `source` молча не срабатывает |
| `env.lac`, `env.w7x` | файлы есть, правил в `get_platform` нет |
| `ASTRA_COMPILER` на любой площадке, кроме TOK | `env.aug_gcc`, `env.darwin_gcc` и т.п. не существуют. Кроме того, `Makefile`/`Makexpr` вызывают `get_platform` с тем же окружением, и на Mac `PLATFORM` станет `darwin_gcc`, т.е. проверка на OpenBLAS не сработает |
| `modules_install.sh` | добавляет суффикс `_${comp}` сам, а `get_platform` — ещё раз, если задана `ASTRA_COMPILER`: возможно имя `tok_gcc_gcc` |
| Linux-ветка `modules_inst_loc.sh` | условие `platform == "linux"` недостижимо (1.2) |
| Версии в `modules_install.sh` и `astra_rc` | установщик кладёт модули в `.../unstable`, а `astra_rc` ищет `nov24`, `may26`, `jun25`, `dec24`, `spider/latest`. После `modules_install.sh` каталоги надо переименовывать вручную или делать симлинки (предположительно так и делается на площадках) |
| `modules_install.sh`, QuaLiKiZ при `comp=gcc` | собирается `TOOLCHAIN=gcc`, но копируется `libQLK-intel-release-default.a` (предположительно ошибка) |
| `modules_install.sh`, ветка gfortran | срабатывает только при `FC` ровно `gfortran` |

### 4.2. Ошибки и опечатки в конкретных файлах

| файл | что | последствие |
|---|---|---|
| `env.cz` | `export PYTHON_BINDIR==/compass/home/tardini/miniforge3/bin` | значение начинается с `=`, каталог в `PATH` недействителен |
| `env.cz` | `INTEL_LICENSE_FILE` из 13 путей с повторами и путём macOS (`/Users/Shared/...`) | безвредно, явный копипаст |
| `env.freia` | `export MPIFC=mpifort`, затем `export MPIFC=mpicc` | `MPIFC=mpicc`, `MPICC` не задан: `Makexpr` (NEO, QuaLiKiZ) компилирует Fortran C-обёрткой, линкует пустым `$(MPICC)` |
| `env.freia`, `env.puhti` | нет `COMP_LIBDIR` | `-L` без пути в `Makexpr` (`-lifcore`) и в `IRC_LIB` для TORBEAM |
| `env.hgw` | `COMP_LIBDIR=${PKG}/compiler/2025.1/bin` (должно быть `lib`) | `LD_LIBRARY_PATH` и `-lirc`/`-lifcore` смотрят не туда (предположительно работает за счёт системных путей) |
| `env.docker` | `export PATH=${PATH}:${COMP_BINDIR}` до определения `COMP_BINDIR` | в `PATH` добавляется пустой элемент (= текущий каталог) |
| `env.docker` | `ASTRA_EXT` = путь AUG `/shares/departments/AUG/...`, `CC="gcc -O"` при `FC=ifort` | вероятно, образ монтирует общий диск AUG |
| `env.sevilla` | `NETCDF="-Wl,-rpath,..."` кладётся в `LD_LIBRARY_PATH` | флаг линковщика в списке каталогов не работает |
| `env.tohtori` | пробел перед `#!/bin/bash`; `COMP_LIBDIR` присвоен дважды (2023.2.0 → 2023.2.1); `COMP_BINDIR` из Intel 2019.4 при библиотеках 2023.2 | компилятор и runtime разных поколений |
| `env.iter` | `echo $LD_LIBRARY_PATH` | печать при каждом `source` (дважды за запуск) |
| `env.w7x`, `env.lac` | нет `FC_FLAGS`/`LD_FLAGS` | сборка без `-r8`, без OpenMP (`-qopenmp`), без `-fPIC` |
| `env.w7x` | `MPIFC=mpif90` при `FC=ifx`; `MPI_COMMAND=" /usr/bin/mpirun -n "` | смешение системного MPI и Intel MPI из `$ASTRA_EXT` |
| `env.lac` | нет `MPICC`, есть только `PYTHON_BIN` | `as_exe` смотрит только на `PYTHON_BINDIR` и возьмёт первый `python3` из `PATH`, а не venv |
| `env.darwin` | в `PATH` добавляются неопределённые `COMP_BIN`, `I_MPI_ROOT/bin` (= `/bin`), `MPI_BIN`; `LD_LIBRARY_PATH` на macOS не действует (нужен `DYLD_*`) | безвредный мусор |

### 4.3. `SL_QUEUE`: два несовместимых формата

С коммита `055b5a46` (2026-02-04) `exe/Build` вставляет `${SL_QUEUE}` в `batch.sh` как есть, а раньше писал `#SBATCH --partition=${SL_QUEUE}`. Состояние:

- `tok_gcc`, `tok_ifx` — новый формат (многострочный блок `#SBATCH ...` и `export OMP_*`);
- `env.tok` — сначала новый блок, а **в конце файла второе присваивание** в старом формате: `SL_QUEUE="p.tok\n#SBATCH --qos=p.tok.48h"`. Оно добавлено в `c23a6d08` (2026-07-22, слияние стеллараторной ветки) и перекрывает первое. В `batch.sh` попадает голая строка `p.tok`. `sbatch` перестаёт читать директивы на первой не-комментарийной строке, поэтому `--qos` игнорируется, а `p.tok` выполняется как команда и падает. Предположительно, на TOK batch-режим сейчас запускается с партицией по умолчанию;
- `env.mit` (`sched_mit_psfc_r8`), `env.tohtori` (`gen05_skylake`) — старый формат, та же проблема.

### 4.4. Жёсткие пути и побочные эффекты

- 11 конфигов указывают `ASTRA_EXT`, Python или netCDF в домашние каталоги конкретных людей (`hgw`, `w7x`, `iter`, `gway`, `cz`, `rat2`, `tohtori`, `freia`, `ncku`, `sevilla`, `lac`; пользователи `tardini`, `g2gtardi`, `tarding`, `papo`, `akadam`, `cooseman`, `alfrodgon`, `astrainstaller`). Остальные — в общие каталоги групп и проектов. Конфиг работает только у того, у кого есть права на чтение.
- `env.lac` при каждом `source` создаёт `$AWD/dat/tcvsimdat/`, читает `$AWD/dirresults` и **делает `rm -r` каталога результатов**, имя которого там записано, после чего создаёт его заново. Файл подключается и в `as_exe`, и в `Build`, т.е. минимум дважды за запуск, а также в `install.sh` и `compReg.sh`. Если `dirresults` нет, `dirname` пустой, и выполняется `mkdir` без аргумента (ошибка, но не удаление).
- `env.gway` экспортирует `PYTHONHOME` глобально — это влияет на любой Python в той же оболочке.
- `env.omega` берёт Intel MPI из установки коммерческого пакета CST Studio Suite 2022 (`/fusion/usc/opt/cst/...`), а компилятор — из Intel 2020.
- `env.perlmutter` использует собственный Intel oneAPI 2022.1 в `/opt/intel/oneapi`, а не Cray PE и не NVIDIA HPC SDK; MPI (2021.10) новее компилятора.
- `env.tok` смешивает Intel MPI (`I_MPI_ROOT`, `mpiifx`) для сборки и OpenMPI 5.0.8 (`MPI_COMMAND`) для запуска VMEC. `NETCDF_LIB` записан в одинарных кавычках, `${NETCDF_HOME}` раскроется только в make.
- `env.tok` использует Python waterboa 2025.06, `tok_ifx`/`tok_gcc` — 2024.06; `tok` и `tok_ifx` в остальном почти совпадают (обе с ifx), но у `tok_ifx` нет `STELLOPT_HOME`, `NETCDF_*`, `MPI_COMMAND`.
- `astra_rc` зашивает путь FFTW MPCDF (`/mpcdf/soft/SLE_15/packages/znver4/fftw/...`) для TORIC на всех площадках.
- Большинство Intel-конфигов **перезаписывают** `LD_LIBRARY_PATH`, а не дополняют его.
- Компиляторы ifort/icc (11 конфигов) Intel убрал из oneAPI 2025; `env.sevilla` гасит именно предупреждение об этом (`-diag-disable=10448`). Эти площадки, видимо, сидят на oneAPI ≤ 2024.
- `Build` перед запуском делает `ulimit -s 100000` (или `unlimited`) — стек нужен большой независимо от площадки.

### 4.5. Полезные флаги

| флаг | где | зачем |
|---|---|---|
| `-r8` | Intel | REAL по умолчанию = 8 байт |
| `-no-prec-div` | Intel | быстрое (менее точное) деление; влияет на побитовое совпадение с GNU-сборкой |
| `-traceback` | Intel | стек вызовов при падении |
| `-qopenmp` / `-fopenmp` | все, кроме `w7x`, `lac` | OpenMP; число потоков IPC ограничивает `MAX_NWORKERS` из `Build` |
| `-diag-disable=10448` | `sevilla` | убрать предупреждение ifort о снятии с поддержки |
| `-ffree-line-length-none`, `-Wl,--no-as-needed` | нет upstream, есть в нашем `env.ubuntu` | длинные строки в `src/astra/a2vmec.f90`; `-lX11` при `--as-needed` (Ubuntu) |
| `I_MPI_FABRICS=shm` | большинство Intel | MPI только внутри узла |
| `OMP_PLACES=cores`, `OMP_PROC_BIND=close` | `tok*` (в `SL_QUEUE`) | привязка потоков на узле 128 ядер |
