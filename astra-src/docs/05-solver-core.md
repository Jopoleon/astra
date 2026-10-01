# 05. Ядро решателя ASTRA 8 (`src/astra`, `src/misc`)

> Состояние исходников: MPCDF `67281112` (2026-10-01).
> Все пути ниже даны относительно корня `astra-src/`. Номера строк — для этого коммита.

## О чём глава

Глава описывает Fortran-ядро ASTRA 8: точку входа `program astra`, инициализацию,
главный цикл по времени, управление шагом, решаемые уравнения переноса (NE, TE, TI,
FP/ток, UPAR, F0..F9), разностную схему и трёхдиагональный решатель, сетку и нормировки,
связь с решателями равновесия через `IPEQL`, хранение переменных, вывод в JSON,
обработку ошибок. Отдельно разобрано, какие места документации `docu/*.tex`
относятся к старым версиям и в коде A8 уже не существуют.

Смежные главы (здесь только ссылки, без подробностей):

- [01-overview.md](01-overview.md) — общая картина и история;
- [02-build-and-run.md](02-build-and-run.md) — сборка, `astra8`, `as_exe.py`, окружение;
- [03-model-language.md](03-model-language.md) — язык модели `equ/<model>` и `pyparse`;
- [04-experiment-data.md](04-experiment-data.md) — формат `exp/<exp>`, U-файлы, `read_input.f90`;
- [06-equilibrium.md](06-equilibrium.md) — FEQIS, SPIDER, EMEQ, `gs_solver`, VMEC подробно;
- [07-heating-and-external-modules.md](07-heating-and-external-modules.md) — NBI, ECRH, TGLF, QLK, RABBIT, TORBEAM, STRAHL;
- [08-graphics-and-output.md](08-graphics-and-output.md) — X11 GUI, JSON→CDF, пост-обработка;
- [09-user-area-catalog.md](09-user-area-catalog.md) — каталог `repo_user_area` (`sbr`, `fnc`, `equ`, `exp`);
- [glossary.md](glossary.md) — переменные, единицы, сокращения.

## Содержание

1. [Архитектура: что такое «ядро» и как оно собирается](#1-архитектура-что-такое-ядро-и-как-оно-собирается)
2. [Файлы `src/astra`, `src/misc`, `inc`](#2-файлы-srcastra-srcmisc-inc)
3. [Хранение переменных: плоские массивы и указатели](#3-хранение-переменных-плоские-массивы-и-указатели)
4. [Точка входа и инициализация](#4-точка-входа-и-инициализация)
5. [Главный цикл по времени и `STEPUP`](#5-главный-цикл-по-времени-и-stepup)
6. [Управление шагом по времени](#6-управление-шагом-по-времени)
7. [Решаемые уравнения](#7-решаемые-уравнения)
8. [Граничные и начальные условия](#8-граничные-и-начальные-условия)
9. [Сетка и нормировки](#9-сетка-и-нормировки)
10. [Численная схема: `RUNEQ`, `TRIDIAG`, `RUNEQ_TETI`](#10-численная-схема-runeq-tridiag-runeq_teti)
11. [Связь с равновесием: `METRIC` и `IPEQL`](#11-связь-с-равновесием-metric-и-ipeql)
12. [Свободная граница и уравнения цепей](#12-свободная-граница-и-уравнения-цепей)
13. [Подпрограммы пользователя внутри шага](#13-подпрограммы-пользователя-внутри-шага)
14. [Вывод: JSON, рестарт, IMAS-подобная структура](#14-вывод-json-рестарт-imas-подобная-структура)
15. [OpenMP и параллельность](#15-openmp-и-параллельность)
16. [Ошибки, остановы и сообщения](#16-ошибки-остановы-и-сообщения)
17. [Расхождения документации `docu/` и кода A8](#17-расхождения-документации-docu-и-кода-a8)
18. [Граф вызовов](#18-граф-вызовов)
19. [Где смотреть в исходниках](#19-где-смотреть-в-исходниках)

---

## 1. Архитектура: что такое «ядро» и как оно собирается

ASTRA — это не одна программа, которая читает модель. **Исполняемый файл
пересобирается под каждую модель.** Модель `equ/<model>` транслируется Python-парсером
`pyparse/parser_main.py` в Fortran, и только потом линкуется вместе с ядром.

Правило линковки в `repo_user_area/exe/Makefile` (цель `$(ATASK)`):

```make
$(PY_TMP_GEN)            # = $(AWD)/pyparse/parser_main.py -equ $(EQU) -exp $(EXP)
$(FC) $(LD_FLAGS) -o $(ATASK) $(ASTRA_SRC) $(TMP_SRC) ... $(LIB) $(X11_FLAG)
```

где

- `ASTRA_SRC = $(ASTRA_DIR)/astra_main.F90 $(WORKFLOW)` — главная программа и «workflow»;
  по умолчанию `WORKFLOW` = `src/astra/stepup.F90` (опция `-workflow` в
  `repo_user_area/exe/as_exe.py` позволяет подставить свой файл вместо штатного `stepup`).
  Поэтому `astra_main.F90` и `stepup.F90` **не входят** в `libastra.a`
  (их нет в списке `AMODS`);
- `TMP_SRC` — сгенерированные парсером файлы в `$(AWD)/src/tmp/`:
  `associate_pointers.f90`, `inivar.f90`, `detvar.f90`, `postep.f90`, `ininam.f90`,
  `init_converge_step.f90`, `eqns_inc.f90` (плюс `declar.fml`, `declar.fnc`,
  `model.txt`, а с X11 — `astra_out.f90`, `set_graph_names.f90`);
- `LIB = -lsbr -lnbi $(X11_LIB) $(QLKNN_LIB) $(TGLF_LIB) $(TRAVIS_LIB) $(TORBEAM_LIB) $(RABBIT_LIB) $(GRAPH_FLAG) -lastra -lfnc -lfeqis -lmisc $(SPIDER_LIB) $(JSON_LIB) -L$(LIB_DIR) $(MKL_LIB) $(RT_LIB) $(TORIC_LIB)`.

Библиотеки ядра (`Makefile`, переменные `MMODS`, `AMODS`, `AOBJS`, `FMODS`, …):

| Библиотека | Источник | Содержимое |
|---|---|---|
| `libmisc.a` | `src/misc/` | `char_manip`, `dbl2char`, `imas_ids`, `numerical_tools`, `pi_const` |
| `libastra.a` | `src/astra/` | все модули из `AMODS` + `ipc_control.o` (C) |
| `libfeqis.a` | `src/feqis/` | решатель равновесия FEQIS (глава 06) |
| `libgraph.a` | `src/graph/` | X11-графика (глава 08) |
| `libnbi.a` | `src/nbi/` | NBI-пакет (глава 07) |
| `libsbr.a`, `libfnc.a` | `$(AWD)/sbr`, `$(AWD)/fnc` | пользовательские подпрограммы и функции (глава 09) |

Порядок модулей в `AMODS` — это и порядок зависимостей при компиляции:
`json_vars → debugger → machine_config → scrunch2d → stella_module → vmec_indata →
parse_utils → emeq → scalars → read_input → time_functions → status → transport_solver →
standard_functions → cpu_usage → parameters_a2equil → surface_contours → flux_avg_imp →
json_rw → gs_solver → a2vmec → a2dkes → metrics → auxiliary → set_x_data`.

Итог: **ядро задаёт каркас** (данные, цикл, шаг, решатель трёхдиагональных систем,
равновесие, I/O), а **физика конкретной модели** (какие уравнения решаются, какие
коэффициенты, какие подпрограммы вызываются) живёт в сгенерированном `eqns_inc.f90` и
соседях. Как пишется модель — глава 03.

## 2. Файлы `src/astra`, `src/misc`, `inc`

В `src/astra` 28 файлов (27 Fortran + 1 C). Ниже все с назначением.

### 2.1. Каркас, данные, цикл

| Файл | Строк | Назначение |
|---|---:|---|
| `astra_main.F90` | 134 | `program astra`: инициализация, начальная сходимость, цикл по времени, нормальный выход |
| `stepup.F90` | 253 | `STEPUP()` — один шаг по времени (workflow по умолчанию) |
| `auxiliary.f90` | 291 | `IFSTEP` (контроль шага), `IFTREQ` (сходимость транспорт↔равновесие), `OLDNEW`, `LINEAV`, `dfj` |
| `scalars.f90` | 222 | модуль скалярных указателей, `scalars_init` (значения по умолчанию) |
| `status.f90` | 459 | модуль профилей (`profiles`, `profiles_x`), `status_init`, `DEFARR`, `SETVAR`, `error_catch` |
| `set_x_data.f90` | 646 | интерполяция эксп. данных в момент `TIME`: `set_x_arrays`, `set_x_scalars`, `get_coil`, `astra_assignments` |
| `json_vars.f90` | 79 | чтение `astra_variables.json` (`read_metadata`, `get_subdict`) |
| `json_rw.f90` | 754 | запись/чтение выходных JSON (`write_ajson`, `read_ajson`, `read_equil`) |
| `cpu_usage.f90` | 91 | счётчики CPU/wall-time, `CPU_report` (печатает `>>> ASTRA normal exit >>>`) |
| `debugger.f90` | 43 | `markloc` — трассировка мест выполнения при `debug>0` |
| `time_functions.f90` | 34 | идентификация вызовов функций со «временной памятью» (`function_id`, `nloc=2200`) |
| `machine_config.f90` | 32 | чтение `exp/cnf/<machine>_description_in.json` (сетка R,Z и т.п.) |
| `parse_utils.f90` | 197 | `assign_val` (разбор `equ/log/<model>`), `path_split`, `split2array`, `inquire_fname` |

### 2.2. Ввод данных

| Файл | Строк | Назначение |
|---|---:|---|
| `read_input.f90` | 1344 | `readInput`, `read_nml` (`tmp/astra.nml`), `read_exp` (разбор `exp/<exp>`), U-файлы (`ufheader`, `ufrd`, `CHECKU`), катушки (`read_coilx`). Подробно — глава 04 |

### 2.3. Перенос

| Файл | Строк | Назначение |
|---|---:|---|
| `transport_solver.f90` | 457 | `RUNEQ` (общее 1D-уравнение), `calc_fxi`, `TRIDIAG`, `GETPEI`, `RUNEQ_TETI`, `TRIDIAG_TETI` |
| `standard_functions.f90` | 1283 | встроенные функции языка модели: `VINT`, `IINT`, `NODE`, `RFA`, `TIMDER`, `TIMAVG`, `FRAMP`, `N_95_POS` и др. |
| `flux_avg_imp.f90` | 300 | `lfs2fsa_impDV` — пересчёт LFS-коэффициентов переноса примеси в усреднённые по поверхности (Angioni 2014, D. Fajardo, дек. 2025) |

### 2.4. Равновесие и геометрия

| Файл | Строк | Назначение |
|---|---:|---|
| `metrics.f90` | 1951 | `METRIC` (диспетчер по `IPEQL`), `EQCYL`, `EQGUESS`, `A2EMEQ`, `A2GSSOLVER`, `BNDRY`, `RHSEQ`, `RHSEQ2`, `CUOFP`, `CUOFMU`, `FPMUOFCU`, `new_grid`, `SETGEO`, `ROC3A`, внешняя метрика |
| `emeq.f90` | 1082 | `EMEQ_MOD` — 3-моментный решатель равновесия (Zakharov), `NP=101` |
| `gs_solver.F90` | 839 | `gssolver` (фиксированная граница, namelist `equil_settings`), `A_EQUIL` (FEQIS/SPIDER) |
| `parameters_a2equil.f90` | 18 | глобальная `equil_now` (`type_equilibrium`) и настройки GS |
| `surface_contours.f90` | 345 | интерполяция контуров поверхностей, `ctr2rz*` (ψ и F на прямоугольной сетке) |
| `scrunch2d.f90` | 254 | `surf2surf`, `scrunch` — Фурье-аппроксимация контура границы |
| `a2vmec.f90` | 529 | интерфейс к VMEC (стеллараторы): запуск, Boozer, `VMECbin2astra` |
| `vmec_indata.f90` | 296 | генерация namelist `&INDATA` для VMEC |
| `a2dkes.f90` | 2343 | интерфейс к DKES-коэффициентам для стеллараторного неоклассического переноса |
| `stella_module.f90` | 9 | общие переменные стеллараторной ветки |

### 2.5. Межпроцессное взаимодействие

| Файл | Строк | Назначение |
|---|---:|---|
| `ipc_control.c` | 252 | семафоры и разделяемая память System V: `initialise_ipc_`, `unlock_sbp_`, `wait4all_`, `fill_int_shm_`, `fill_dbl_shm_`, `sbp2astra_`. Используется обёртками `repo_user_area/sbr/{tglf,qlk,neo}_ipc.f90` для параллельного запуска TGLF/QuaLiKiZ/NEO как дочерних процессов |

### 2.6. `src/misc` и `inc`

| Файл | Назначение |
|---|---|
| `src/misc/pi_const.f90` | `GP=3.14159265359`, `GP2=2GP`, `GP2_sq`, `muvac=4e-7·GP`, `mu0=0.4·GP` (последний — в «единицах ASTRA», MA и м) |
| `src/misc/imas_ids.f90` | IMAS-подобные типы `type_equilibrium` (`profiles_1d`, `eqgeometry`, `coord_sys`, `global_param`) и `equil_allocate` (в коммите `67281112` добавлено выделение `rho_tor_norm`) |
| `src/misc/numerical_tools.f90` | интерполяция (`linterp`, `qinterp`, `reinterp_back`, сплайны), `SMOOTH`, `DERIV`, `INTEGR`, `GRADIENT`, `EXTRAP`, сортировка, `mom2rz` |
| `src/misc/char_manip.f90` | строки: `clean_string`, `to_upper/lower`, `split_string`, `str_in_list` |
| `src/misc/dbl2char.f90` | `ISNUM`, `fmt_smart`, `itoa` |
| `inc/Astra.h` | заголовок C для `ipc_control.c`: `sys/ipc.h`, `sys/sem.h`, `sys/shm.h` и `union semun` (не на Apple) |

## 3. Хранение переменных: плоские массивы и указатели

### 3.1. Идея

Все именованные переменные ASTRA (`TE`, `NE`, `BTOR`, `CF1`, …) — это Fortran-**указатели**
на элементы нескольких плоских массивов. Порядок и имена задаются файлом
`astra_variables.json` в корне. Это позволяет:

- писать/читать JSON целиком по группам (`write_ajson` / `read_ajson`), не перечисляя переменные;
- генерировать код парсером, зная индекс каждой переменной.

### 3.2. Массивы

| Группа JSON | Массив | Модуль | Размер | Примеры |
|---|---|---|---:|---|
| `variables` | `varValues` | `scalars` | 135 | `AB`, `BTOR`, `IPL`, `RTOR`, `ZRD1..96`, `GVAC`, `PHIEDG` |
| `variables_x` | `varxValues` | `scalars` | 135 | те же с суффиксом `X` (значения из эксп. файла) |
| `constants` | `constValues` | `scalars` | 162 | `CF1..16`, `CV1..16`, `CHE1..4`, `CNB1..4`, …, `CDHJ9` |
| `control` | `controlValues` | `scalars` | 15 (+4·`n_sbr`) | `DPOUT`, `TIME`, `TAUMIN`, `TAUMAX`, `TAUINC`, `DELVAR`, `TEND`; хвост — `DTEQ(4,n_sbr)` |
| `internInt` | `internIntValues` | `scalars` | 33 | `NEQUIL`, `INUME1..4`, `IPART`, `IPEQL`, `NA1`, `NA1N..NA19` |
| `internDbl` | `internDblValues` | `scalars` | 46 | `TAU`, `HRO`, `ROC`, `FTO`, `PSIFB`, `ATREQ`, `RBDOT`, `BBDOT` |
| `profiles` | `profiles(NRD, n_prof)` | `status` | 454 | `TE`, `NE`, `FP`, `G11`, `CAR1..128`, `F0..F9` |
| `profiles_x` | `profiles_x(NRD, n_profx)` | `status` | 106 | `TEX`, `NEX`, `CAR1X..64X`, `G11X` |

`NRD = 801` — максимальное число радиальных узлов (`status.f90`).

Связывание делает сгенерированный `src/tmp/associate_pointers.f90`, например
(в рабочей установке `~/astra/a8/src/tmp/associate_pointers.f90`):

```fortran
CF1  => constValues(1)
DTEQL => controlValues(12)
DTEQ(1:4, 1:n_sbr) => controlValues(n_control+1: n_control + 4*n_sbr)
ZEF  => profiles(:, 441)
FTPT => profiles(:, 454)
```

Там же выставляется число пользовательских подпрограмм `n_sbr` и выделяются `DTEQ`, `TEQ`.

Переменные, **не входящие** в JSON (служебные целые/массивы): `KEY`, `NSTEPS`, `ITREQ`,
`LEQ(NEQNS=19)`, `NITOT`, `TEQ(:)`, `DTEQ(4,n_sbr)`, `exp_header`, `tbeg_eq/tend_eq`
(`scalars.f90`).

### 3.3. Чтение метаданных

`json_vars::read_metadata` (`src/astra/json_vars.f90:48`) открывает файл
`'astra_variables.json'` **в текущем каталоге** (то есть в `$AWD` при запуске) через
json-fortran и строит списки имён по группам `get_subdict`. Имена усечены до
`character(len=6)` — длинных имён переменных в ядре быть не может.
Ошибка разбора печатается (`'Error: ', error_msg`), но явного останова в этом месте нет.

Тот же файл читают `pyparse/equ_parser.py`, `pyparse/exp_parser.py`, `pyparse/code_gen.py`
(генерация `associate_pointers`). Структура файла описана в [glossary.md](glossary.md#astra_variablesjson).

### 3.4. Значения по умолчанию

`scalars_init` (`src/astra/scalars.f90:101`) обнуляет все группы, кроме
`constValues = 1.` (все `CF*`, `CV*`, `CHE*`, … по умолчанию равны 1), и задаёт:

| Параметр | Значение | Параметр | Значение |
|---|---|---|---|
| `DPOUT` | 0.01 с | `IPEQL` | 5 (FEQIS) |
| `TAUMIN` | 1e-6 с | `INUME1..4` | 22 |
| `TAUMAX` | 0.05 с | `ITFBE` | 1e6 (FBE фактически выключено) |
| `TAUINC` | 1.1 | `NA1`=`NB1`=`NAB` | 41 |
| `DELVAR` | 0.1 | `ATREQ` / `PTREQ` | 1e-4 / 1e-5 |
| `TEND`, `TPAUSE` | 1e5 с | `AB`=`ABC` | 0.3 м |
| `BTOR` | 3 Т | `RTOR` | 1.5 м |
| `IPL` | 0.3 МА | `AMJ` / `ZMJ` | 2 / 1 |
| `GN2E`, `GN2I` | 0 | `WNE`=`WTE`=`WTI` | 0.03 м |
| `GVAC` | 4.5 | `PHIEDG` | +∞ (`ieee_value`, маркер «не задано») |
| `LEQ(:)` | −1 | `TEQ(:)` | −1e3 |
| `DTEQ(1..4,:)` | 0, −99999, 99999, −1 | `TAU` | 1e-6 |

Профили по умолчанию задаёт `status_init` (`src/astra/status.f90:90`): `TE=TI=0.01`,
`NE=NI=0.1`, `CU=0.1`, `MU=0.3`, `FP` линейный и т.п. Это лишь «страховка» до
`INIVAR`.

## 4. Точка входа и инициализация

`program astra` (`src/astra/astra_main.F90`). Последовательность вызовов:

| Строка | Вызов | Что делает |
|---:|---|---|
| 42–43 | `CPU_TIME`, `SYSTEM_CLOCK` | старт таймеров |
| 45 | `read_metadata` | читает `astra_variables.json` |
| 46 | `associate_pointers` | (сгенерирован) указатели → плоские массивы |
| 47 | `cpu_init` | таймеры подпрограмм |
| 49 | `scalars_init` | значения по умолчанию |
| 50 | `ininam` | (сгенерирован) `DTEQ`, имена подпрограмм `sbr_name`, `LEQ(j)`, подписи графиков из модели |
| 52 | `readInput` | namelist `tmp/astra.nml`, `equ/log/<model>`, `exp/<exp>` (глава 04) |
| 53 | `allocate(CCOIL, VCOIL)` | токи/напряжения катушек |
| 56–60 | `set_graph_names`, `gui_init` | только с X11 и если `TASK(1:3) /= 'BGD'` |

Дальше две ветки.

**Рестарт** (`restart > 0`, строки 63–74): `read_ajson(restart)` загружает все группы из
ранее записанного JSON, `set_x_arrays(2)`, для `IPEQL==5` — `transport2fbe_init`
(передача параметров в FEQIS), `TEND`/`TPAUSE` берутся из namelist.

**Обычный старт** (строки 75–107):

```text
status_init                ! страховочные профили
astra_assignments          ! геометрия, сетки, NA1, ROC, начальная метрика (цилиндроподобная)
IPART = 1                  ! «секция начальных итераций»
set_x_arrays(1)            ! эксп. профили в TSTART без интерполяции по времени
INIVAR                     ! (сгенерирован) NE=NEX, TE=TEX, ..., MU, CC
SETVAR                     ! ZEF=max(1,ZEFX), NI=NE/ZMJ, проверка ABC<=AB
DETVAR                     ! (сгенерирован) присваивания из модели (IPEQL, NEQUIL, TAUMAX...)
EQGUESS                    ! аналитическая метрика по форме границы
INIVAR
transport2fbe_init         ! если IPEQL==5
do while (jt_req == 0)     ! цикл согласования профилей и равновесия
    set_x_scalars; DETVAR; DEFARR; set_x_arrays(1); INIVAR
    NITOT = NITOT + 1
    INIT_CONVERGE_STEP     ! (сгенерирован) присваивания модели без решения уравнений
    IFBEY = 0              ! FBE в инициализации запрещено
    METRIC                 ! равновесие
    jt_req = IFTREQ(ATREQ)
enddo
```

Цикл инициализации продолжается, пока `IFTREQ` не вернёт ненулевое значение.
При `IPART == 1` в `IFTREQ` число итераций ограничено `ITREQMAX = 200`
(`src/astra/auxiliary.f90:143,156`); при превышении печатается
`>>> Warning >>> Max number of initialization iterations is overwhelmed, stopping` и
выполняется `stop` (строки 169–173). Критерий сходимости — максимум относительного
изменения `G11`, `G22`, `VR`, `FP` по узлам и `ROC` между итерациями, сравниваемый с
`ATREQ` (строки 176–212). В сгенерированном `detvar.f90` модели часто переопределяют
`ATREQ` (в тестовой установке `ATREQ = 0.5`).

`astra_assignments` (`src/astra/set_x_data.f90:488`) — важное место:

- проверяет `NA1 <= NRD` (`>>> FATAL ERROR: The radial grid size out of range...`, `ERROR STOP`);
- выставляет `IFDFVX = 4` («запрещено менять») для `AB`, `AWALL`, `RTOR`, `ELONM`, `TRICH`;
- `ABC`/`ELONG` из заданной границы, `ROC = ROC3A(...)`, `FTO`, `FTN`, `BTN`, `ROWALL`;
- строит сетки (`XRHO`, `SXHO`, `RHO`, `SRHO`, раздел 9), `SETGEO`, `NEW_GRID`;
- начальная «цилиндроподобная» метрика `G22`, `VR`, `VRS`, `G11`, `VOLUM`;
- `PSIAX`, `PSIBO`, `FP_NORM`, `rho_pol`, копии `*O`, `UPS0`, `TIMEQL`, `TAUPRP = TAUMIN`;
- `>>> Error: ABC cannot exceed AB` → `ERROR STOP`; предупреждения про `AWALL`.

## 5. Главный цикл по времени и `STEPUP`

### 5.1. Цикл в `astra_main.F90` (строки 122–132)

```fortran
t_prev_json = TSTART - DPOUT
do while (TIME < TEND + 1.d-8)
    if ((TIME - t_prev_json + 1.d-8) > DPOUT) then
        call write_ajson()
        t_prev_json = t_prev_json + DPOUT
    endif
    call STEPUP()
enddo
call CPU_report('>>> ASTRA normal exit >>>')
STOP
```

Следствия:

- первый JSON пишется в `TSTART`, далее примерно каждые `DPOUT` секунд модельного
  времени; запись происходит **до** шага, поэтому состояние после последнего шага в
  JSON может не попасть (если `TEND` не совпадает с сеткой `DPOUT`) — *предположительно*,
  специально не проверялось;
- `TPAUSE` в пакетном цикле не используется (встречается только в GUI-коде,
  *предположительно*; в `astra_main.F90` только присваивается при рестарте);
- `IFSTEP` ограничивает шаг сверху значением `DPOUT` (раздел 6), так что моменты записи
  не «перепрыгиваются» больше чем на шаг.

### 5.2. Один шаг: `STEPUP` (`src/astra/stepup.F90`)

| Строки | Действие |
|---:|---|
| 29–30 | `BTN = BTOR`, `FTN = FTO` (значения на начало шага) |
| 34–42 | на первом шаге `tau_old = tau_new = tau`; `tau = tau_old`; `IPART = 2` («эволюция») |
| 46 | `NIO = NI` |
| 49–56 | переключение на свободную границу: если `TIME >= ITFBE`, то `IFBEY += 1`, иначе `IFBEY = 0` |
| 58 | `detvar()` |
| 63–66 | если `plasma_up .or. ifbey==0`: `OLDNEW` (сохранить `*O`), `set_x_scalars` |
| 73–77 | сброс `PSIFBO`, `RBDOT`, `BBDOT`, `PSIEXO`, `PSPLXO` |
| 79–85 | токи катушек: из уравнений цепей (`CCOIL`) при `TIME > ITFBE`, иначе из эксп. трасс (`get_coil`) |
| 87–95 | логика `IBCPSI` (ψ как граничное условие, до 3) |
| 97 | **`time_step_accuracy: do`** — повтор шага с меньшим `tau`, если точность не достигнута |
| 99–109 | `tau = tau_new`, `ITREQ = 0`, сброс граничных величин |
| 113 | **`tr_eq_loop: do while (jkey == 0)`** — цикл «перенос ↔ равновесие» |
| 118–120 | `set_x_arrays(2)` — эксп. профили с интерполяцией по времени и на новую метрику |
| 122 | `METRIC()` — равновесие/метрика (раздел 11) |
| 124–133 | `RBDOT = (FTO−FTN)/(FTO+FTN)/TAU`, `BBDOT = (BTOR−BTN)/(BTOR+BTN)/TAU`; обнуляются при `IFBEY == 1` |
| 135–138 | `hro is nan` → `error_catch` |
| 140–142 | **`eqns_inc(IBCPSI, bc_type_for_fp, dfpdrbm12)`** — решение уравнений модели (только если `plasma_up`) |
| 145–176 | проверки `te/ti/ne/fp isinf` (среднее \|x\| > 1e8) и `isnan` → `error_catch` |
| 179–183 | GUI-клавиши; `jkey = IFTREQ(ATREQ)` |
| 186–199 | `PSIFB = PSIEXT − PSPLEX·ROC·dfpdrbm12`; сдвиг `FP` при `ITFBP<0`, `IBCPSI==1` |
| 203–205 | `DEFARR()` — проверки положительности, продление профилей за `NA1` |
| 207–215 | `tau_old = tau`; если `IFSTEP()==0` и `IFBEY /= 1` — повтор с новым `tau`, иначе выход |
| 219–223 | `tau` зажимается в `[TAUMIN, TAUMAX]`, `tau_new = tau`, `tau = tau_old` |
| 226–236 | обновление цепей после шага: `SPIDUPDATE` (`IPEQL==4`, только `#ifdef SPIDER`) или `FEQISUPDATE` (`IPEQL==5`), если `IFBEY>=1` и `ICIRCQ>0` |
| 238–240 | `TIME = TIME + TAU`, `TAUPRP = tau`, `POSTEP()` (сгенерирован, подпрограммы с `>`) |
| 242–246 | если `POSTEP` изменил `TAU` (например `TSCTRL`), он учитывается как следующий шаг |

**Важно:** во время эволюции (`IPART = 2`) `IFTREQ` возвращает 2 сразу: в коде
`NTREQ = 1` жёстко (`src/astra/auxiliary.f90:154–158`). Значит, `tr_eq_loop`
выполняется **ровно один раз** на попытку шага — итераций перенос↔равновесие внутри шага
нет, метрика берётся «с запаздыванием» на шаг (комментарий в `IFTREQ`: «iterations
switched off (timing VR(t) downshifted)»). Комментарий в `stepup.F90:115`
(«NITREQ regulates this») устарел: переменной `NITREQ` в коде нет.

## 6. Управление шагом по времени

`IFSTEP` (`src/astra/auxiliary.f90:22–113`):

$$
C_\tau = \max\left(\frac{1}{\mathtt{TAUINC}},\ \max_j \frac{|y^{old}_j/y_j - 1|}{\mathtt{DELVAR}}\right),
$$

где максимум берётся по `NE`, `TE`, `TI` (если `LEQ(1..3) > 0`) и `F0..F9` (если
`LEQ(10..19) > 0`; для `Fj` используется `dfj`). Затем

$$
\tau_{new} = \max\left(\mathtt{TAUMIN},\ \min\left(\mathtt{TAUMAX},\ \tau/C_\tau,\ \mathtt{DPOUT}\right)\right).
$$

- если `τ_new >= τ` — шаг принят (`IFSTEP = 1`, `NSTEPS += 1`), `τ_new` станет следующим шагом;
- иначе шаг отвергнут: все `NE, NI, TE, TI, FP, VR, UPAR, F0..F9` восстанавливаются из `*O`,
  `ROC = ROCO`, `FTO = FTN`, вызывается `CUOFP`, `IFSTEP = 0` — шаг повторяется с `τ_new`.

Смысл: шаг растёт не быстрее чем в `TAUINC` раз за шаг, а относительное изменение
любой эволюционирующей величины за шаг не должно превышать `DELVAR`. `FP`, `UPAR`
в критерии **не участвуют**. Шаг с `IFBEY == 1` (первое включение свободной границы)
принимается безусловно (`stepup.F90:209`).

Легенда `LEQ(j)` (комментарий в `auxiliary.f90:25–34`): `-1` — уравнения нет,
`0` — `AS` (присваивание), `1` — `EQ` (уравнение), `2` — `FU`, `3` — поток тепла +
конвекция. Индексы: 1 `NE`, 2 `TE`, 3 `TI`, 4 `CU`, 5 равновесие, 10..19 — `F0..F9`.
Значения `LEQ` выставляет сгенерированный `ininam.f90`.

В моделях часто задают `TAUMIN = TAUMAX` (фиксированный шаг; так в тестовой
установке, `detvar.f90`: `TAUMIN=TAUMAX=5e-2`).

## 7. Решаемые уравнения

Источники: `docu/section3.tex` (общая формулировка), `docu/section4.tex`
§«Transport equations in Astra notations» (строки 173–380), `docu/mom_transp_eq_add_EF.tex`
(момент), `docu/stiff_transport_EF.tex` (жёсткий перенос), генератор
`pyparse/eqns.py`, `pyparse/const_text.py`.

### 7.1. Операторы

С $\rho = \sqrt{\Phi/(\pi B_0)}$, $V' = \partial V/\partial\rho$ (section4.tex:179–193):

$$
\hat D_n[f] = \frac{1}{V'}\left(\partial_t - \frac{\dot B_0}{2B_0}\partial_\rho\rho\right)[V' f],\qquad
\hat D_T[f] = \frac{3}{2}\frac{1}{(V')^{5/3}}\left(\partial_t - \frac{\dot B_0}{2B_0}\partial_\rho\rho\right)\left[(V')^{5/3} f\right],
$$

$$
\hat D_\psi[f] = \partial_t f - \frac{\rho\dot B_0}{2B_0}\partial_\rho f .
$$

Член с $\dot B_0$ — адиабатическое сжатие; в коде он представлен величинами `BBDOT`
и `RBDOT` (раздел 10).

### 7.2. Основная система (уравнение `Amain`, section4.tex:195–227)

$$
\begin{aligned}
\hat D_n[\mathtt{NE}] + \partial_V(\mathtt{QN}) &= \mathtt{SN} + \mathtt{SNN}\cdot\mathtt{NE},\\
\hat D_T[\mathtt{NE}\cdot\mathtt{TE}] + \partial_V\left(625\,\mathtt{QE} + \tfrac52\gamma_e\,\mathtt{TE}\cdot\mathtt{QN}\right) &= 625\,(\mathtt{PE} + \mathtt{PET}\cdot\mathtt{TE}),\\
\hat D_T[\mathtt{NI}\cdot\mathtt{TI}] + \partial_V\left(625\,\mathtt{QI} + \tfrac52\gamma_i\tfrac{\mathtt{NI}\,\mathtt{TI}}{\mathtt{NE}}\mathtt{QN}\right) &= 625\,(\mathtt{PI} + \mathtt{PIT}\cdot\mathtt{TI}),\\
2\pi\rho\,\mathtt{CC}\,\hat D_\psi[\mathtt{FP}] + V'(\mathtt{CUBS} + \mathtt{CD}) &= 5\,\mathtt{IPOL}^2\,\partial_\rho\left[\mathtt{G22}\,\partial_\rho \mathtt{FP}\right].
\end{aligned}
$$

- $625 = 1/(1.6\cdot10^{-3})$ — из-за смешанных единиц (кэВ·10¹⁹ м⁻³ против МДж);
- $\mathtt{GN2E} = \tfrac52\gamma_e$, $\mathtt{GN2I} = \tfrac52\gamma_i$, по умолчанию 0 (и в доке, и в `scalars_init`);
- `PE`, `PI` — явные источники, `PET`, `PIT` — коэффициенты при `TE`, `TI` (неявная часть).

### 7.3. Матрица переноса (section4.tex:258–300)

$$
\begin{pmatrix}
\dfrac{\mathtt{QN} - \mathtt{SLAT}\cdot\mathtt{GNX}}{\mathtt{G11}\cdot\mathtt{NE}}\\[2mm]
\dfrac{625\,\mathtt{QE}}{\mathtt{G11}\cdot\mathtt{NE}\cdot\mathtt{TE}}\\[2mm]
\dfrac{625\,\mathtt{QI}}{\mathtt{G11}\cdot\mathtt{NI}\cdot\mathtt{TI}}\\[2mm]
\dfrac{0.4\pi}{B_p}\mathtt{CUBS}
\end{pmatrix}
= -
\begin{pmatrix}
\mathtt{DN} & \mathtt{HN} & \mathtt{XN} & \mathtt{CN}\\
\mathtt{DE} & \mathtt{HE} & \mathtt{XE} & \mathtt{CE}\\
\mathtt{DI} & \mathtt{HI} & \mathtt{XI} & \mathtt{CI}\\
\mathtt{DC} & \mathtt{HC} & \mathtt{XC} & 0
\end{pmatrix}
\begin{pmatrix}
\partial_\rho \ln \mathtt{NE}\\ \partial_\rho \ln \mathtt{TE}\\ \partial_\rho \ln \mathtt{TI}\\ -1
\end{pmatrix}.
$$

Диагональные элементы — `DN` (диффузия частиц), `HE` (χₑ), `XI` (χᵢ); `CN`, `CE`, `CI` —
конвективные скорости; `DC`, `HC`, `XC` — бутстреп-коэффициенты. Генератор `pyparse/eqns.py`
(`neeqn`, `teeqn`, `tieqn`) раскладывает недиагональные члены в явные/неявные слагаемые
коэффициентов `A`, `B` для `RUNEQ`.

### 7.4. Метрика (section3.tex:348–410)

$$
\mathtt{G11} = V'\langle(\nabla\rho)^2\rangle,\quad
\mathtt{G22} = \frac{R_0}{J}\frac{V'}{4\pi^2}\left\langle\frac{(\nabla\rho)^2}{r^2}\right\rangle,\quad
\mathtt{G33} = \left\langle\frac{R_0^2}{r^2}\right\rangle,
$$

$\mu = \mathtt{MU} = 1/q = \frac{1}{2\pi B_0\rho}\partial_\rho\psi$,
ток $I_{pl} = \frac{G_2}{\mu_0}\partial_\rho\psi$; в коде
`IPL = 5·IPOL(NA1)·G22(NA)·(FP(NA1)−FP(NA))/HRO/(2π)/RTOR`
(`pyparse/eqns.py`, `cueqn`). `IPOL` — нормированный полоидальный ток $J = I/(R_0B_0)$.

### 7.5. Уравнение тока / полоидального потока

Генерируется `pyparse/eqns.py::cueqn` + шаблоны `pyparse/const_text.py::CUEQN`. В
`eqns_inc.f90` (рабочая установка) видно:

- проводимость `CC` (в тестовой модели — формула Sauter с `FTLLMR`), бутстреп `CUBS`,
  внешний ток `CD` (например `CUBM + CUECR`);
- `YWD = (CUBS + CD)·YD/(IPOL³·G33)`, `YWB = 0.4π·CC·YC/IPOL²`;
- граничное условие выбирается по тому, что определено в модели:
  - задан `IPL` (без `UEXT`/`LEXT`) → `bctype = 3`, `bc_values(5) = HRO·μ0/G22(NA)·IPL·RTOR/IPOL(NA1)`
    (для `IPEQL ∈ {6,7,9}` — стеллараторная формула с `SG11`, `MV`);
  - задан `UEXT` → `FP(NA1) = FPO(NA1) + TAU·UEXT`, `bctype = 1`;
  - задан `LEXT` → уравнение цепи: `PSIEXT`, `PSPLEX = LEXT/ROC·5·IPOL(NA1)·G22(NA)/(2π)/RTOR`, смешанное условие;
- после решения: `UPL = (FP−FPO)/TAU − YQDCMF`, `ULON = IPOL·G33·(UPL − 2π·ROC²·BTOR·BBDOT·MU)`, `CUOFP`
  (пересчёт `CU`, `MU`), `PSIAX`, `PSIBO`, `FP_NORM`, `rho_pol`.

### 7.6. Вспомогательные уравнения `F0..F9` (section4.tex `Aux`, `QFj`)

$$
\frac{1}{V'}\partial_t(V'\,\mathtt{Fj}) + \partial_V(\mathtt{QFj}) = \mathtt{SFj} + \mathtt{SFFj}\cdot\mathtt{Fj},\qquad
\mathtt{QFj} = \mathtt{G11}\left(\mathtt{VFj}\cdot\mathtt{Fj} - \mathtt{DFj}\,\partial_\rho\mathtt{Fj}\right).
$$

В A8 их 10 (`F0..F9`, индекс `LEQ(10..19)`), в старой документации — 3. Генератор —
`pyparse/eqns.py::fjeqn` (вызов `RUNEQ` со `SFFj`, `SFj`). Часто используются для примесей,
быстрых частиц, моделей пьедестала.

### 7.7. Тороидальный момент `UPAR`

`docu/mom_transp_eq_add_EF.tex` (Astra 7.0, строки 100–175): переменная
$\mathtt{UPAR} = \langle U_\parallel B\rangle/B_0$, уравнение

$$
\frac{1}{V'}\partial_t\big(V'\langle M_\phi\rangle\big) + \frac{1}{V'}\partial_\rho\left[V' R_\phi - V'\rho\frac{\dot B_0}{2B_0}\langle M_\phi\rangle\right] = \left\langle R\sum_s S_{\phi,s}\right\rangle .
$$

В коде (`pyparse/eqns.py::upeqn`, `pyparse/const_text.py::UPEQN`): коэффициент при
производной по времени `UPS0 = MRHO·G41·RTOR/IPOL`, неоклассические поправки `UPS1`
(`DLNEO`), `UPS2` (`SGNEO`), диффузия `XUPAR`, пинч `CNPAR`, остаточные напряжения
`RUPAR/RUPYR/RUPFR`, источник момента `TTRQ` (явный) и `TTRQI`, численный флаг `INUME4`,
граничный узел `NA1U`. Поправки `UPS1/UPS2` включаются при `IPROT >= 1`.
Замечание: в `upeqn` присваивается `YWC(J) = 0.0*TTRQ(J)`, но затем в `UPEQN.eqn`
`TTRQ` передаётся в `RUNEQ` как явный источник, а `TTRQI` — как `YWD`; *предположительно*,
`YWC` здесь не используется (не проверено построчно в сгенерированном коде).

### 7.8. Нейтралы

В документации (section3.tex:1393–1445) — кинетическое уравнение в приближении
пластины; в section4.tex:2621 — подпрограмма `NEUT`. В A8 ядро нейтралов **не решает**:
`NEUT` — пользовательская подпрограмма (`repo_user_area/sbr/`), вызываемая из
`eqns_inc.f90` по расписанию `DTEQ` (раздел 13); в тестовой модели `sbr_name(1)="NEUT"`,
`DTEQ(1,1) = 5e-2`. Её результаты (`NN`, `TN`, `SNNBM` …) входят в источники как обычные
профили. Подробнее — глава 07/09.

### 7.9. Жёсткий перенос (`DVN`, `DVE`, `DVI`, `DSE`…)

`docu/stiff_transport_EF.tex` (Pereverzev & Corrigan, CPC 2008): к диффузии добавляется
искусственная `DV`, а к конвекции — компенсирующий член, так что на сошедшемся решении
поток не меняется, но схема устойчива при резкой зависимости χ от градиента. Реализовано
**в генераторе**, а не в `RUNEQ` (заголовок `transport_solver.f90:13` прямо говорит:
«Grigori's scheme with additional D and V for stiff transport is not yet implemented»):
`pyparse/eqns.py` (≈ строки 263–267, 389–450, 669–742) добавляет `DV` к `A` и член
`YWB = DV·log(y_j/y_{j+1})/HRO` к `B`. Артефакты этого видны в профилях
`SDN`, `PDE`, `PDI`, `SD0..SD9` («Artificial source introduced by DV*»).

## 8. Граничные и начальные условия

### 8.1. На оси

Нулевой поток: в `RUNEQ` коэффициент `CC(1) = 0` (`transport_solver.f90:119`), то есть
уравнение в первом узле не связано с несуществующим узлом 0.

### 8.2. На краю

Типы `bctype` в `RUNEQ` (`transport_solver.f90:52–55`):

| `bctype` | Смысл | Задаётся в модели |
|---:|---|---|
| 1 | значение на границе `y(Ngridb)` | `NEB`, `TEB`, `TIB`, `FjB`, или «внешнее условие» `NE = ...` за `ND1` |
| 2 | поток через границу `Q(b)` | `QNB`, `QEB`, `QIB`, `QFjB` |
| 3 | смешанное `C5·y(b) + C6·y(b−1) = C7` | ток `IPL`, уравнение цепи `LEXT`, ψ от FBE |
| 4 | поток пропорционален значению `Q = c·y` | в `pyparse/` не генерируется (grep `bctype = 4` пуст); поддержан только в `RUNEQ` |

Логика выбора — `pyparse/eqns.py::bnd_init`:

- граница расчёта `ND1` = `NA1` или `NODE(ρ_BC)` если в модели `XX:EQ[ρ]` / `EQ[j,ρ,tbeg,tend]`
  (функция `rhoBC`); за `ND1` профиль задаётся выражением модели или `XXX` (эксп.),
  а при отсутствии внешнего условия — линейной интерполяцией к сепаратрисе;
- одновременно задать значение и поток нельзя (`Too many boundary conditions`);
- `QNB` может задаваться как `QNNB·NE` (`QNNB`, `QETB`, `QITB`, `QFFjB` в `internDbl`);
  в шапке `RUNEQ` предупреждение: «at the moment Qb is explicit, no option for QNNB, QETB, QITB».

Документация: section3.tex:1107–1228, section4.tex:555–783 (Table 4.8).

### 8.3. Начальные условия

Section3.tex:1021–1105 описывает три варианта для тока: задать `j₀`, `μ` или
стационарное `U_pl = const`. В A8 начальные условия задаёт сгенерированный `inivar.f90`
(например `NE = NEX`, `TE = TEX`, `TI = TIX`, `UPAR` из `VTORX`, `MU = 1/(CAR2X + 1e-4)`),
после чего цикл инициализации (раздел 4) согласует их с равновесием. `INIT_CONVERGE_STEP`
выполняет присваивания модели без решения уравнений переноса.

## 9. Сетка и нормировки

Источник: `docu/norm_grid_add_EF.tex`, код `set_x_data::astra_assignments`,
`metrics::new_grid` (`metrics.f90:1799`).

| Величина | Формула | Смысл |
|---|---|---|
| `x` | $\sqrt{\Phi/\Phi_b}$ | нормированная координата |
| `ρ` | $\sqrt{\Phi/(\pi B_0)}$, $B_0 = $ `BTOR` | «размерная» координата, м |
| `ROC` | $\sqrt{\mathtt{FTO}/(\pi\,\mathtt{BTOR})}$ | ρ на границе (`metrics.f90:44–46`) |
| `HROX` | $1/(\mathtt{NA1} - 0.5)$ | шаг по `x` |
| `XRHO(j)` | $(j - 0.5)\cdot\mathtt{HROX}$ | основная сетка (центры) |
| `SXHO(j)` | $j\cdot\mathtt{HROX}$ | сдвинутая сетка (грани) |
| `RHO`, `SRHO`, `HRO` | `XRHO·ROC`, `SXHO·ROC`, `HROX·ROC` | размерные аналоги |

Таким образом, первый узел не на оси, а на `HROX/2`, а последний основной узел ровно
на границе: `XRHO(NA1) = (NA1−0.5)·HROX = 1`; узел `NA1` отстоит от `NA` на полный шаг, а от
грани `SXHO(NA)` — на половину.
Потоки (`QN`, `QE`, …) и коэффициенты `A`, `B`, `G11` живут на сдвинутой сетке.
Комментарий в `TRIDIAG` (`transport_solver.f90:225–226`): «boundary value is assumed to be
on the main last grid point, so 1−dx/2. This has to be corrected later on...».

`NA1 ≤ NRD = 801`; `NB1` — до стенки (`AWALL`), `NAB` — до `AB`. Профили продлеваются
за `NA1` экспоненциально с ширинами `WNE`, `WTE`, `WTI` (`DEFARR`).

При изменении `ROC` (геометрия меняется во времени) сетка по `ρ` «дышит», и в уравнениях
появляются члены адиабатического сжатия с
`RBDOT = (FTO−FTN)/(FTO+FTN)/TAU` ≈ $\dot\Phi_b/(2\Phi_b)$ и
`BBDOT = (BTOR−BTN)/(BTOR+BTN)/TAU` ≈ $\dot B_0/(2B_0)$ (`stepup.F90:125–126`).

## 10. Численная схема: `RUNEQ`, `TRIDIAG`, `RUNEQ_TETI`

### 10.1. Общее уравнение `RUNEQ`

Заголовок `src/astra/transport_solver.f90:33–42`:

$$
\frac{1}{G}\partial_t(G H y) + \frac{1}{V}\partial_x\Big(c\,(Q - \ldots)\Big) = S\,y + P,\qquad
Q = -G_{11}\left(A\,\partial_x y + B\,y\right) + G_{11} R,
$$

с членами адиабатического сжатия, добавляемыми к `B` и `S` через `rbdot`, `bbdot`.
`G`, `V`, источники — на основной сетке; `A`, `B`, `R`, `G11` — на сдвинутой.
`R` ненулевой только в уравнении момента (остаточное напряжение).

Сигнатура (строки 8–9):

```fortran
RUNEQ(G_new, H_new, G_old, H_old, y_old, W_in, V_in, unit_coeff, G11, A_in, B_in, R_in,
      Src_y_in, Src_in, rbdot, bbdot, Ngridb, Ngrid, dx, dt, x_in, imethod, bctype,
      bc_values, y_out, Q_out, adcmp_term, mphit)
```

### 10.2. Выбор схемы `imethod` (= `INUME1..4`)

| `INUMEn` | Модуль | θ | `f(ξ)` (`calc_fxi`, строки 165–200) |
|---:|---|---:|---|
| 21 | неявная, центральные разности | 1 | $1 + \xi/2$ |
| 22 | неявная, степенной закон (power law) | 1 | $(1 - 0.1|\xi|)^5_+ + \max(\xi,0)$ |
| 23 | неявная, экспоненциальная | 1 | $\xi/(1 - e^{-\xi})$ |
| 31 | Кранк–Николсон, центральные | 0.5 | как 21 |
| 32 | Кранк–Николсон, степенной | 0.5 | как 22 |
| 33 | Кранк–Николсон, экспоненциальная | 0.5 | как 23 |

$\xi_j = \Delta x\,B_j/\max(A_j, 10^{-16})$ — сеточное число Пекле; $g(\xi) = f(\xi) - \xi$.
Номера: `INUME1` — `NE`, `INUME2` — `TE`/`TI`, `INUME3` — `FP`, `INUME4` — `UPAR`.
По умолчанию 22.

Явные схемы 11/12, упомянутые в описании `INUME1` в `astra_variables.json`, **в
`RUNEQ` не реализованы**: `SELECT CASE(imethod)` знает только 21–23 и 31–33, и при другом
значении `theta` остаётся неинициализированной (строки 94–99).

### 10.3. Дискретизация

Для внутренних узлов (строки 110–128) собирается трёхдиагональная система
$A_j y_{j+1} + B_j y_j + C_j y_{j-1} = R_j$:

$$
\begin{aligned}
A_j &= -\theta\,g_v\,v_j\,f_j, \qquad C_j = -g_{v,j}\,v_{j-1}\,g_{j-1},\\
B_j &= (GH)^{new}_j + \theta\left(-\Delta t\,G\,S_y + g_{v,j}(v_j g_j + v_{j-1}f_{j-1})\right),\\
R_j &= (GH)^{old}_j y^{old}_j + \Delta t\,G\,P + (1-\theta)\left(\ldots y^{old}\ldots\right),
\end{aligned}
$$

где $g_v = \Delta t\,G/(\Delta x^2 V)$, $v = c\,G_{11}A$. Граничное условие на краю
подставляется по `bctype` (строки 130–143), затем `TRIDIAG` (прогонка, строки 203–237).
Поток `Q_out` восстанавливается из решения на сдвинутой сетке; в граничной точке —
экстраполяцией (`bctype` 1, 3), значением (`2`) или `c·y` (`4`).

Источник от сжатия (строки 84–91):
`adcmp_term = bbdot·Pdot_1 + (rbdot − bbdot)·Pdot_2` — добавляется к явному источнику и
возвращается наружу (используется для `FP`, переменная `YQDCMF`).

### 10.4. Неявный обмен энергией: `RUNEQ_TETI`

Если в модели уравнения записаны как `TE*: EQ` / `TI*: EQ` (stiff_transport_EF.tex),
генератор `pyparse/eqns.py::tetieqn` вызывает совместный решатель `RUNEQ_TETI`
(`transport_solver.f90:276`) с блочной прогонкой `TRIDIAG_TETI` (строка 385): обмен
энергией $P_{ei}(T_e - T_i)$ входит неявно в обе системы.

Коэффициент обмена `GETPEI(j)` (строка 240):

- если модель задаёт `PEIQI` ≠ 0 — используется он;
- иначе кулоновская формула:
  $P_{ei} = 0.00246\,\ln\Lambda\,n_e n_i Z^2/(A\,T_e^{3/2})$, $\ln\Lambda = 15.9 - 0.5\ln n_e + \ln T_e$;
  при `|IPROT| ∈ {2, 4}` суммирование идёт по основным ионам и трём примесям (`SUZPEI`).
- побочный эффект: `GETPEI` правит `PET`, `PIT` для корректного вывода `PETOT`/`PITOT`.

### 10.5. Ньютон и итерации

`docu/newton_solver.tex` описывает ньютоновский решатель (`IPNWT`, файл
`dat/newton_numerics.dat`). **В коде A8 его нет**: ни `IPNWT`, ни `newton_numerics` не
встречаются в `src/`, `pyparse/`, `repo_user_area/exe/`. Нелинейность обрабатывается
только (а) уменьшением шага `IFSTEP` и (б) схемой жёсткого переноса (раздел 7.9).

## 11. Связь с равновесием: `METRIC` и `IPEQL`

`METRIC` (`src/astra/metrics.f90:16–253`). При `IPART ∈ {1,3}` пересчитывает
`ROC = sqrt(FTO/π/BTOR)`. Дальше диспетчер по `IPEQL` (строки 118–228):

| `IPEQL` | Действие | Подпрограммы |
|---:|---|---|
| −2 | цилиндр, без тороидальности | `EQCYL`, `RHSEQ` |
| −1 | метрика из эксп. файла (`G11X`, `G22X`, `VRX` …) | `ROC3A`, `set_external_metric`, `RHSEQ` |
| −6 | метрика из внешнего файла (не X-массивы) | `set_external_metric_2`, `RHSEQ` |
| 0 | без решателя: аналитическая оценка по форме границы | `EQGUESS` |
| 1 | 3-моментный решатель EMEQ | `RHSEQ`, `A2EMEQ` → `EMEQ_MOD` |
| 3 | только правые части (p′, FF′, `CUTOR`) | `RHSEQ` |
| 4 | SPIDER (закрытый код, `#ifdef SPIDER`) | `RHSEQ2`, `A2GSSOLVER(3)` |
| 5 | **FEQIS** (по умолчанию) | `RHSEQ2`, `A2GSSOLVER(101)` |
| 6 | стелларатор, метрика из внешнего файла | `set_external_metric_2`, `RHSEQ` |
| 7 | стелларатор, VMEC | `EQCYL_STELLA`, `RHSEQ`, `vmec_interface` |
| 2, 9, прочие | ветки нет — метрика не меняется | — |

Замечания:

- **Интервал пересчёта** — `DTEQL`: решатель вызывается, если `TIME − TIMEQL ≥ DTEQL`
  (на первом шаге — всегда). `DTEQL = 0` (по умолчанию) — каждый шаг.
- `IPEQL = 9` упоминается в стеллараторных ветках `CUOFP` и генераторе (`IPEQL == 6 .or. 9 .or. 7`),
  но в `METRIC` ветки для 9 нет.
- Для `IPEQL < 3` после расчёта метрики заполняется IMAS-подобная `equil_now`
  (`psi`, `F_dia`, контуры `R,Z` по формуле с `SHIF`, `AMETR`, `TRIA`, `ELON`), строки 228–251, —
  чтобы JSON-вывод имел равновесные группы и в этих режимах.
- Ошибка EMEQ на начальных итерациях: `Equilibrium problem at the initial iterations` → `ERROR STOP`
  (строки 149–156).
- Время, проведённое в равновесии, копится в `wallTime_equ` и печатается в `CPU_report`
  строкой `Equilibrium`.

`A2GSSOLVER` (`metrics.f90:1099`) собирает вход для `gs_solver::gssolver`:
профили `FP`, давление (с вкладом пучка при `NB2EQL`), граница из `BNDRY`
(эксп. группы `BND`/`BNDU`, 8-точечная 3-моментная граница или параметрическая), напряжения
катушек `VCOIL`, вращение `omega_rot = VTOR/(RTOR+SHIF+AMETR)`. Результат — новая метрика
`G11, G22, G33, VR, VRS, IPOL, SLAT, …`, форма (`SHIF`, `ELON`, `TRIA`, `AMETR`) и новый `ROC`;
затем `new_grid`. `gssolver` читает namelist `equil_settings` из `nml_file`
(`exp/nml/<exp>` или `exp/nml/<machine>`) и описание машины
`exp/cnf/<machine>_description_in.json`. Подробности решателей — глава 06.

Подозрительное место (*не проверено в рантайме*): в ветке с `GSSOLVER`
(`metrics.f90:1222–1292`) локальная переменная `equil_out` не заполняется, но после
вызова из неё берутся `PSPLEX = equil_out%global_param%psplex` и
`PSIEXT = −2π·equil_out%global_param%psiext`. Поля имеют инициализацию по умолчанию
`-9.0D40` (`src/misc/imas_ids.f90`, `type_global_param`), так что, *предположительно*,
`PSPLEX`/`PSIEXT` получают мусорные значения. На фиксированной границе они используются
только в `PSIFB` (`stepup.F90:186–189`), которое применяется лишь при `ITFBP < 0`.

## 12. Свободная граница и уравнения цепей

Управляющие величины (`astra_variables.json`, `stepup.F90`):

| Имя | Смысл |
|---|---|
| `ITFBE` | момент включения FBE (free boundary equilibrium); отрицательное — выключено; по умолчанию 1e6 |
| `IFBEY` | счётчик шагов с FBE (0 — фиксированная граница/PBE) |
| `ITFBP` | 0 — ничего; ≠0 — использовать `PSIFB` как граничное условие для ψ; `<0` — неявно |
| `IBCPSI` | счётчик режима «ψ как граничное условие» (до 3) |
| `ICIRCQ` | 1 — решать уравнения цепей катушек (`FEQISUPDATE` / `SPIDUPDATE` после шага) |
| `IPCTRL` | контроллер внешних катушек |
| `CCOIL`, `VCOIL` | токи/напряжения катушек: до `ITFBE` — из эксп. групп `CCOILX`/`VCOILX` (`get_coil`), после — из цепей |
| `PSIFB`, `PSIEXT`, `PSPLEX` | ψ на границе, внешний поток, «Green-функционный» поток от плазмы |
| `IPLFBE`, `PTREQ` | ток из FBE и допуск по току |
| `plasma_up` | логический флаг в `metrics` (по умолчанию `.true.`); `.false.` — пробой, транспорт не решается |

При первом шаге с FBE (`IFBEY == 1`) `RBDOT`, `BBDOT` обнуляются (геометрия меняется
скачком), а шаг принимается без проверки точности. Если `plasma_up = .false.` и FBE
включено, `A2GSSOLVER` вызывает `feqis_main` напрямую с током `IPLFBE`
(`metrics.f90:1213–1220`).

Известная остановка FEQIS: `Error in find new boundary (totpoints > nz2*nr2)` в
`src/feqis/fbe_core.f90:705` с `stop` (код возврата 0). На регрессии `fbe/AUG33040_2500`
она происходит около t = 2.62 (по монорепо `AGENTS.md`). Разбор — глава 06.

## 13. Подпрограммы пользователя внутри шага

Пользовательские подпрограммы из модели (глава 03) вставляются генератором в разные места:

- **внутри `eqns_inc`** — по расписанию. Шаблон (рабочая установка `eqns_inc.f90:51–60`):

  ```fortran
  IFSUB = 0
  if (KEY /= 0 .and. ABS(KEY - DTEQ(4,1)) < 0.1) IFSUB = 1
  if (TIME >= DTEQ(2,1) .and. TIME <= DTEQ(3,1) .and. TIME - TEQ(1) + 1.E-7 - DTEQ(1,1) > 0.0) IFSUB = 1
  if (IFSUB == 1) then
      TEQ(1) = TIME
      call NEUT(...)        ! с учётом времени в wallTime_sbr(1)
  ```

  `DTEQ(1,j)` — интервал вызова, `(2,j)`/`(3,j)` — окно по времени, `(4,j)` — клавиша GUI,
  `TEQ(j)` — время последнего вызова;
- **«<»** — до решения уравнений (`detvar`, в начале шага);
- **«>»** — после шага, в `POSTEP` (`stepup.F90:240, 248–251`). Если такая подпрограмма меняет
  шаг (например `TSCTRL`), это учитывается в строках 242–246.

Время каждой подпрограммы печатается в отчёте `CPU_report` строками `Subroutine <name>`.

## 14. Вывод: JSON, рестарт, IMAS-подобная структура

`write_ajson` (`src/astra/json_rw.f90:387`):

- файл `$AWD/ncdf_out/<exp><equ>-<n>.json`, где `n = j_call + restart` (счётчик вызовов);
- группы: `variables`, `variables_x`, `constants`, `control`, `internInt`, `internDbl`,
  `profiles_x(1:NA1)`, `profiles(1:NA1)`, `equil_signals` (8 скаляров из `equil_now%global_param`),
  `equil_profiles` (28), `equil_rz2d`, `equil_coord`, `equil_rect`;
- имена ключей берутся из тех же списков `astra_variables.json`, что и указатели;
- печать `   Written file ...`.

`read_ajson(n_restart)` (строка 722) и `read_equil` (строка 235) читают такой файл
обратно — это механизм рестарта (`restart` в `tmp/astra.nml`, раздел 4).

Конвертация JSON → CDF делается **вне Fortran** функцией `json2cdf` в
`repo_user_area/exe/as_exe.py` (глава 08).

`equil_now` (`src/astra/parameters_a2equil.f90`) — глобальная структура
`type_equilibrium` из `src/misc/imas_ids.f90`: это не настоящий IMAS IDS, а упрощённая
копия схемы ITM/IMAS (`profiles_1d`, `eqgeometry%boundary`, `eqgeometry%rectgrid`,
`coord_sys`, `global_param`). Значения-«пустышки» — `-9.0D40`. Служит обменным форматом
между ASTRA и решателями равновесия и источником равновесных групп JSON.

## 15. OpenMP и параллельность

- В `src/` **нет** ни одной директивы `!$omp` и нет `use omp_lib`. Распараллеливание
  ядра через OpenMP отсутствует.
- Флаги OpenMP есть только в окружении сборки: `platform/env.tok_gcc`, `platform/env.darwin`
  (`FC_FLAGS="-O0 -w -fPIC -fopenmp"`, `LD_FLAGS="-fopenmp"`), `-qopenmp` в Intel-окружениях,
  `OMP_FFLAG` в `repo_user_area/exe/Makexpr`. В рабочей установке `~/astra/a8/platform/env.ubuntu`
  (создаётся инсталлятором монорепо, в `astra-src` его нет) —
  `FC_FLAGS="-O0 -w -fPIC -fopenmp -ffree-line-length-none"`, `LD_FLAGS="-fopenmp -Wl,--no-as-needed"`.
  Флаги нужны внешним модулям и BLAS/LAPACK.
- Реальная параллельность: (а) многопоточный BLAS (MKL/OpenBLAS; в монорепо рекомендовано
  `OPENBLAS_NUM_THREADS=1`), (б) дочерние процессы TGLF/QuaLiKiZ/NEO через System V IPC
  (`src/astra/ipc_control.c`): `initialise_ipc_` создаёт набор семафоров и сегменты
  разделяемой памяти, пишет файл `tmp/<exp><equ>.ipc`, запускает подпроцессы;
  `unlock_sbp_`/`wait4all_` синхронизируют; `fill_*_shm_`/`sbp2astra_` передают данные;
  (в) MPI внутри VMEC (`MPI_COMMAND` в `a2vmec.f90`).

## 16. Ошибки, остановы и сообщения

**Коды возврата.** Большинство «ошибок времени счёта» делают `stop` (код 0), а не
`ERROR STOP`: `error_catch`, `IFTREQ` (лимит итераций), `DEFARR` при `MU <= 0`, `ROC3A`,
FEQIS `find new boundary`. Вместе с тем, что `as_exe.py` всегда возвращает 0
(монорепо `AGENTS.md`), **успех нужно определять по строке `>>> ASTRA normal exit >>>`**
и наличию выходных файлов, а не по коду возврата.

| Сообщение | Где | Тип |
|---|---|---|
| `>>> ASTRA normal exit >>>` + таблица времени | `astra_main.F90:131` → `cpu_usage.f90:26` | норм. выход |
| `hro is nan`; `te/ti/ne/fp isinf`; `te/ti/ne/fp isnan` → `TE Fp NE G11 ...` `somethings not right, quit run` | `stepup.F90:135–176`, `status.f90:446` | `stop` |
| `>>> ERROR >>> Time = ...` `The rotational transform is less or equal zero at rho_N = ...` | `status.f90:204–214` (`DEFARR`) | `stop` |
| `>>> ERROR >>> Time = ... DEFVAR` `The variable "NE" is less or equal zero` (и `ZEF`, `NI`, `TE`, `TI`) | `status.f90:224–251` | `ERROR STOP` |
| `>>> Warning >>> Max number of initialization iterations is overwhelmed, stopping` | `auxiliary.f90:169–173` | `stop` |
| `Equilibrium problem at the initial iterations` | `metrics.f90:149–156` | `ERROR STOP` |
| `>>> Error >>> ROC3A >>> Illegal input: R+Delta < a` | `metrics.f90:808–817` | `STOP` |
| `>>> Warning >>> Maximum size of the equilibrium grid is 101` | `metrics.f90:901` (EMEQ, `NP`) | предупр. |
| `>>> FATAL ERROR: The radial grid size out of range. Parameter "NA1" cannot exceed 801` | `set_x_data.f90:511–513` | `ERROR STOP` |
| `>>> Error: ABC cannot exceed AB. Check your file exp/...` | `set_x_data.f90:533–534` | `ERROR STOP` |
| `>>> Warning >>> Inconsistent boundary setting.` | `status.f90:435` (`SETVAR`) | предупр. |
| `>>> read_equ_log: Error, empty model file name`; `>>> Error: file "equ/log/<m>" missing` | `read_input.f90:87–102` | `ERROR STOP` |
| ошибки формата `exp`, U-файлов, `Unknown input type`, `GRIDTYPE ... not implemented` | `read_input.f90` (≈40 мест) | `ERROR STOP` (глава 04) |
| `>>> READAT: File "..." reading error` | `parse_utils.f90:31–33` | `ERROR STOP` |
| `Too many time functions calls` | `time_functions.f90` (`nloc = 2200`) | останов |
| `Error in find new boundary (totpoints > nz2*nr2)` | `src/feqis/fbe_core.f90:705` | `stop` |

Особенность проверки положительности в `DEFARR` (`status.f90:216–224`): сравнивается
**максимум** профиля по радиусу (`YNE = max(YNE, NE(j))`), поэтому ошибка возникает, только
если величина неположительна **во всех** узлах. Локальные отрицательные значения эта
проверка не ловит. Сообщение `The 3M equilibrium solver does not converge` (ветка `'3M'`)
при этой логике недостижимо.

Трассировка: `debugger::markloc` при `debug ≥ 1..3` (параметр `debug` в `tmp/astra.nml`)
печатает метки `Trackback:` с именами мест — основной инструмент локализации падения.

## 17. Расхождения документации `docu/` и кода A8

| Документация | Что сказано | Что в коде A8 |
|---|---|---|
| `newton_solver.tex` | `IPNWT`, `dat/newton_numerics.dat` | отсутствует |
| `norm_grid_add_EF.tex` | `FLXDR`, `ADCMPF` | отсутствуют; сжатие всегда через `RBDOT`/`BBDOT` |
| `norm_grid_add_EF.tex`, `stepup.F90:115` | `NITREQ` — число итераций Tr-Eq | переменной нет; `NTREQ = 1` жёстко в `IFTREQ` |
| `norm_grid_add_EF.tex` | `INUME = 22` — «centered» | 22 = неявная, степенной закон (`calc_fxi`) |
| описание `INUME1` в JSON | 11/12 — явные схемы | не реализованы в `RUNEQ` |
| section4.tex:785–945, описание `NEQUIL` | выбор решателя по `NEQUIL` (≤0 нет, ≤41 3M, >41 ESC) | выбор по `IPEQL`; ESC в A8 нет; `NEQUIL`/`MEQUIL` — размеры сеток |
| section4.tex | `F1..F3` | `F0..F9` |
| язык модели (`EQ0`, `E3M`, `ESC`, `SPIDER`, `EQS` …) | выбор решателя равновесия в модели | `pyparse/equ_parser.py:165–185` пишет код в `LEQ(5)`, но ядро `LEQ(5)` не читает; решает только `IPEQL` (обычно из `detvar`) |
| section3/4 | нейтралы как часть системы | пользовательская `NEUT` |
| описание `IPROT` в JSON | флаг полоидальных членов в моменте | также переключает формулу `GETPEI` (`|IPROT| = 2, 4`) |

## 18. Граф вызовов

```text
program astra                                    src/astra/astra_main.F90
├─ read_metadata                                 json_vars.f90
├─ associate_pointers                            src/tmp (gen)
├─ cpu_init, scalars_init, ininam(gen)
├─ readInput                                     read_input.f90
│   ├─ read_nml (tmp/astra.nml), config_read(machine)
│   ├─ assign_val (equ/log/<model>)              parse_utils.f90
│   └─ read_exp (exp/<exp>, udb/ U-файлы)
├─ [X11] set_graph_names, gui_init
├─ restart>0: read_ajson → set_x_arrays(2) → transport2fbe_init
└─ иначе:
    ├─ status_init, astra_assignments            status.f90, set_x_data.f90
    ├─ set_x_arrays(1), INIVAR, SETVAR, DETVAR, EQGUESS, INIVAR
    ├─ transport2fbe_init (IPEQL==5)
    └─ do while IFTREQ==0  (≤200)
        set_x_scalars, DETVAR, DEFARR, set_x_arrays(1), INIVAR,
        INIT_CONVERGE_STEP(gen), METRIC, IFTREQ
do while TIME < TEND
├─ write_ajson (каждые DPOUT)                    json_rw.f90
└─ STEPUP                                        stepup.F90
    ├─ detvar(gen); OLDNEW; set_x_scalars; get_coil
    └─ time_step_accuracy:
        ├─ tr_eq_loop (1 проход):
        │   ├─ set_x_arrays(2)
        │   ├─ METRIC                            metrics.f90
        │   │   ├─ EQCYL | set_external_metric(_2) | EQGUESS
        │   │   ├─ RHSEQ + A2EMEQ → EMEQ_MOD      emeq.f90
        │   │   ├─ RHSEQ2 + A2GSSOLVER → gssolver → A_EQUIL → feqis_main | spider_run
        │   │   └─ vmec_interface → VMEC, Boozer, VMECbin2astra   a2vmec.f90
        │   ├─ eqns_inc(gen)
        │   │   ├─ sbr по DTEQ (NEUT, NBI, TGLF/QLK через IPC, ...)
        │   │   ├─ NE:  RUNEQ                     transport_solver.f90
        │   │   ├─ TE,TI: RUNEQ | RUNEQ_TETI(GETPEI)
        │   │   ├─ FP:  RUNEQ → CUOFP
        │   │   ├─ UPAR: RUNEQ
        │   │   └─ F0..F9: RUNEQ
        │   ├─ проверки inf/nan → error_catch
        │   └─ IFTREQ (→2)
        ├─ DEFARR
        └─ IFSTEP (повтор при 0)
    ├─ FEQISUPDATE | SPIDUPDATE (FBE + ICIRCQ)
    ├─ TIME += TAU
    └─ POSTEP(gen)
CPU_report('>>> ASTRA normal exit >>>'); STOP
```

## 19. Где смотреть в исходниках

| Файл | Что там |
|---|---|
| `src/astra/astra_main.F90` | точка входа, инициализация, цикл по времени, нормальный выход |
| `src/astra/stepup.F90` | шаг по времени: FBE-логика, Tr-Eq цикл, проверки, контроль шага |
| `src/astra/auxiliary.f90` | `IFSTEP` (шаг), `IFTREQ` (сходимость, `NTREQ=1`, `ITREQMAX=200`), `OLDNEW`, `LINEAV` |
| `src/astra/transport_solver.f90` | `RUNEQ`, схемы `INUME`, `TRIDIAG`, `GETPEI`, `RUNEQ_TETI` |
| `src/astra/scalars.f90` | скалярные указатели, `scalars_init` (все умолчания) |
| `src/astra/status.f90` | профили, `status_init`, `DEFARR` (проверки), `SETVAR`, `error_catch` |
| `src/astra/set_x_data.f90` | интерполяция эксп. данных, `astra_assignments` (сетки, `NA1`) |
| `src/astra/metrics.f90` | `METRIC`/`IPEQL`, EMEQ/GS-интерфейсы, `CUOFP`, `new_grid`, `RHSEQ*` |
| `src/astra/emeq.f90` | 3-моментный решатель равновесия |
| `src/astra/gs_solver.F90` | `gssolver`, `A_EQUIL` → FEQIS/SPIDER |
| `src/astra/a2vmec.f90`, `vmec_indata.f90`, `a2dkes.f90`, `stella_module.f90` | стеллараторная ветка (VMEC, Boozer, DKES) |
| `src/astra/surface_contours.f90`, `scrunch2d.f90` | контуры поверхностей, Фурье-аппроксимация границы |
| `src/astra/flux_avg_imp.f90` | LFS→FSA пересчёт коэффициентов переноса примеси |
| `src/astra/standard_functions.f90` | встроенные функции модели (`VINT`, `IINT`, `TIMDER`, …) |
| `src/astra/read_input.f90`, `parse_utils.f90`, `machine_config.f90` | ввод модели, эксп. данных, конфигурации машины |
| `src/astra/json_vars.f90`, `json_rw.f90` | метаданные и JSON-вывод/рестарт |
| `src/astra/cpu_usage.f90`, `debugger.f90`, `time_functions.f90` | тайминги, трассировка, функции со «временной памятью» |
| `src/astra/parameters_a2equil.f90`, `src/misc/imas_ids.f90` | `equil_now`, IMAS-подобные типы |
| `src/astra/ipc_control.c`, `inc/Astra.h` | System V IPC для TGLF/QLK/NEO |
| `src/misc/numerical_tools.f90`, `pi_const.f90`, `char_manip.f90`, `dbl2char.f90` | численные и строковые утилиты, константы |
| `astra_variables.json` | реестр всех переменных (порядок = индексы массивов) |
| `pyparse/eqns.py`, `pyparse/const_text.py` | генерация уравнений, BC, жёсткий перенос, ток, момент |
| `repo_user_area/exe/Makefile` | сборка библиотек и исполняемого файла под модель |
| `~/astra/a8/src/tmp/*.f90` (рабочая установка) | пример сгенерированных `eqns_inc`, `inivar`, `detvar`, `associate_pointers` |
| `docu/section3.tex`, `docu/section4.tex` | физическая постановка, единицы, переменные |
| `docu/norm_grid_add_EF.tex`, `stiff_transport_EF.tex`, `mom_transp_eq_add_EF.tex`, `newton_solver.tex` | сетка, жёсткий перенос, момент, (отсутствующий) Ньютон |
