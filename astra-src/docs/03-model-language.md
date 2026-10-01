# 03. Язык модели ASTRA 8 и парсер `pyparse`

## О чём глава

Глава описывает файл модели `equ/<model>`: его синтаксис в ASTRA 8 (A8), отличия от ASTRA 7 (A7) и то, как Python-парсер `pyparse/` превращает модель в Fortran-код `src/tmp/*.f90`, который затем компилирует `make`. Разобраны все модули парсера, сгенерированные файлы и их место в цикле счёта. Пример построчно сверен с рабочей установкой `~/astra/a8` (модель `flux_feqis`, эксперимент `aug34954`). В конце собраны известные ловушки парсера.

Состояние исходников: `astra-src`, MPCDF `67281112` (2026-10-01). Пути даны относительно корня `astra-src`. Если вывод получен прогоном парсера в копии под `/tmp`, это сказано явно. Если вывод получен только чтением кода, он помечен «по коду».

Связанные главы: [01-overview.md](01-overview.md), [02-build-and-run.md](02-build-and-run.md) (сборка и запуск), [04-experiment-data.md](04-experiment-data.md) (файлы эксперимента, `X`-массивы), [05-solver-core.md](05-solver-core.md) (решатель переноса, `LEQ`, `ND1`), [06-equilibrium.md](06-equilibrium.md) (`IPEQL`, FEQIS), [07-heating-and-external-modules.md](07-heating-and-external-modules.md) (подпрограммы `sbr/`), [08-graphics-and-output.md](08-graphics-and-output.md) (вывод `\` и `_`), [09-user-area-catalog.md](09-user-area-catalog.md) (каталог `equ/`, `fml/`, `fnc/`, `sbr/`), [glossary.md](glossary.md) (переменные и единицы). Установка: [`../../installer/docs/DETAILED.md`](../../installer/docs/DETAILED.md).

## Содержание

1. [Место парсера в сборке](#1-место-парсера-в-сборке)
2. [Лексика файла модели](#2-лексика-файла-модели)
3. [Классы строк](#3-классы-строк)
4. [Имена: переменные, массивы, формулы, функции](#4-имена-переменные-массивы-формулы-функции)
5. [Выражения и их перевод в Fortran](#5-выражения-и-их-перевод-в-fortran)
6. [Уравнения переноса: `TE:EQ`, `:AS` и граничные условия](#6-уравнения-переноса)
7. [Подпрограммы: `NAME(args)<:dt:t1:t2:key`](#7-подпрограммы)
8. [Вывод: радиальные профили `\` и временные сигналы `_`](#8-вывод)
9. [Отличия A8 от A7](#9-отличия-a8-от-a7)
10. [Модули `pyparse` по отдельности](#10-модули-pyparse-по-отдельности)
11. [Сгенерированные файлы и их место в счёте](#11-сгенерированные-файлы-и-их-место-в-счёте)
12. [Разбор примера: `flux_feqis` → `src/tmp`](#12-разбор-примера-flux_feqis--srctmp)
13. [Ловушки и ошибки парсера](#13-ловушки-и-ошибки-парсера)
14. [Как проверять модель без сборки](#14-как-проверять-модель-без-сборки)
15. [Где смотреть в исходниках](#15-где-смотреть-в-исходниках)

---

## 1. Место парсера в сборке

Порядок запуска подробно описан в [02-build-and-run.md](02-build-and-run.md). Для языка модели важно следующее:

| шаг | кто | что происходит |
|---|---|---|
| 1 | `exe/as_exe` → `repo_user_area/exe/as_exe.py` | пишет `tmp/astra.nml` (namelist `&astra_log`: `exp_file`, `equ_file`, `tbeg_nml`, `tend_nml`, `machine` …), вызывает `exe/Build` |
| 2 | `exe/Build` | `make -f exe/Makefile MODEL=… FEXP=… ATASK=…` |
| 3 | `repo_user_area/exe/Makefile` | правило `PY_TMP_GEN = $(AWD)/pyparse/parser_main.py -equ $(EQU) -exp $(EXP)` пишет `src/tmp/*` |
| 4 | `make` | компилирует `src/tmp/*.f90` вместе с ядром, `libsbr` (из `sbr/`) и `libfnc` (из `fnc/`) в `bin/<model>_batch.exe` / `_gui.exe` |
| 5 | исполняемый файл | читает `tmp/astra.nml`, `equ/log/<model>`, **сам** читает файл эксперимента (`src/astra/read_input.f90::read_exp`) |

Следствия:

- Модель транслируется в Fortran во время сборки. Каждая модель даёт свой исполняемый файл.
- Файл эксперимента парсер читает только для двух целей: выдать предупреждения о `X`-массивах, которых нет в эксперименте, и вычислить `nr_x_max`, максимальное число радиальных точек в профилях эксперимента (`pyparse/parser_main.py`, строки 17–33). Сами данные читает Fortran во время счёта (см. [04-experiment-data.md](04-experiment-data.md)).
- В `Makefile` цель `$(ATASK)` зависит от `$(EQU)`, `$(PYFILES)`, `$(FMLS)` и исходников ядра, но **не** от `$(EXP)` (`repo_user_area/exe/Makefile`). Если поменять только файл эксперимента, `src/tmp` не перегенерируется и `nr_x_max` в `associate_pointers.f90` остаётся старым. Безопасный способ — `touch equ/<model>` или `make clean`. Приводит ли устаревший `nr_x_max` к реальной ошибке — не проверено: прогон `30000_3.4` (сетка 201 точка) бинарником, собранным с `nr_x_max=51`, за 10 с не упал.
- Если в модели есть ошибка, парсер часто только пишет сообщение в лог, а `as_exe` всегда возвращает код 0 ([DETAILED.md §3.9](../../installer/docs/DETAILED.md)). Смотрите вывод.

## 2. Лексика файла модели

Предобработку выполняет `pyparse/parse_as.py::equ_prepare` (строки 118–153):

| правило | реализация |
|---|---|
| начальные и конечные пробелы обрезаются, пустые строки пропускаются | `line.strip()`, `if line == ''` |
| блочный комментарий: строка, **начинающаяся** с `/*`, открывает блок, строка, **кончающаяся** на `*/`, закрывает | строки 140, 149. В docstring функции ошибочно сказано про `%` |
| строка, начинающаяся с `!`, — комментарий | строка 143 |
| текст после первого `!` отбрасывается | `line.split('!')[0]`. Это действует и внутри кавычек |
| `;` разделяет операторы в одной строке | `newline.split(';')`. Это тоже действует внутри кавычек |
| строка, кончающаяся на `@!@`, включает `config.checkeqn` | тогда `eqns_init` добавляет `NI=F1+…+F9` (модуль `iondens_ass`, см. §10) |

Затем `EQU_PARSER` (`pyparse/equ_parser.py`, строки 90–102) раскладывает строки по классам:

- строки с `=` или `:` идут в `equ_lines`, и из них **удаляются все пробелы** (`line.replace(' ','')`);
- строки с `\` идут в радиальный вывод (`radout`) без удаления пробелов;
- строки с `_`, в которых нет `:`, `\` и `=`, идут во временной вывод (`timout`);
- остальные строки (например, `Profile output`) никуда не попадают и фактически являются комментариями. Это правило A7: строка без спецсимволов — комментарий (`docu/section5.tex`, «Setting a model»).

Регистр не важен: строки переводятся в верхний регистр (`tmp.upper()`, строка 123). Поэтому строковые аргументы подпрограмм тоже оказываются в верхнем регистре (§13).

Длина строки парсером не ограничена. Генератор режет Fortran-строки длиннее 110 символов продолжением `&` (`parse_as.py::write_fortran`, строки 16–51). В A7 предел был 132 символа (`section5.tex`).

## 3. Классы строк

| вид строки | пример | обработка |
|---|---|---|
| присваивание переменной или профилю | `CF3=NA1;` `HE=HCKIM;` | `detv_lines` → `detvar.f90` или, для коэффициентов и потоков уравнений, `eqns_lines` → `eqns_inc.f90` |
| тип уравнения | `TE:AS;` `CU:EQ;` `NE:EQ[1,0.8]` | `assign_d[key]` (`equ_parser.py`, строки 141–152) |
| неявная связка TE–TI | `TE*:EQ;` | `assign_d['TE'] = 'implicit_EQ'` → общий решатель `tetieqn` (`eqns.py`, строка 645) |
| подпрограмма | `NEUT:.05:::;` `TORBA<:.05:::;` | `sbr_lines` → `parse_sbr` (§7) |
| радиальный вывод | `Tex\TE\\TEX;` | `radout` → `set_graph_names.f90`, `astra_out.f90` |
| временной вывод | `q95_1/MU(AFX(0.95));` | `timout` → там же |
| вставка Fortran «как есть» | `"WORK(J,3)=TE(J)"` | в радиальный цикл `detvar.f90` без трансляции (`code_gen.py`, ветка `line[0]=='"'`) |
| `:=` и `=>` | `CF1:=TE(AFX(0.))` | заменяются на `=` (строка 123) |

**Повторное присваивание.** Если левая часть уже встречалась, запись хранится под ключом `KEY|n` (строки 127–134). Порядок операторов сохраняется, оба присваивания попадают в код. Так работает, например, `PBOL1 = …; PBOL1 = PBOL1-PBR1` в `flux_feqis`.

**Куда попадает присваивание.** Если левая часть входит в `config.eqn_list` или в списки коэффициентов и потоков (`config.coeff_d`, `config.flux_d`), строка идёт в блок уравнений. Иначе она идёт в `DETVAR` (`equ_parser.py`, строки 136–139). Списки заданы в `pyparse/config.py`:

| уравнение | коэффициенты (`coeff_d`) | потоки/источники (`flux_d`) | граничные условия (`bnd_d`) |
|---|---|---|---|
| `TE` | `DE HE XE CE DVE DSE` | `PE PET` | `TEB QEB QETB` |
| `TI` | `DI HI XI CI DVI DSI` | `PI PIT` | `TIB QIB QITB` |
| `NE` | `DN HN XN CN DVN DSN` | `SN SNN` | `NEB QNB QNNB` |
| `CU` | `DC HC XC CC CD` | `MV MU` | `IPL MUB` (ключ `MU`) |
| `UPAR` | `CNPAR XUPAR RUPAR RUPYR RUPFR` | `TTRQ` | `UPARB TTRQB` |
| `F0`…`F9` | `DFj VFj GFj DVFj DSFj` | `SFj SFFj` | `FjB QFjB QFFjB` |

`config.eqn_list` = `NE TE TI CU Equil UPAR F0 … F9`, всего 16 элементов → `LEQ(1..16)` (`config.py`). Ключ `Equil` записан в смешанном регистре, а строки модели приведены к верхнему, поэтому строка `Equil:…` никогда не распознаётся как уравнение (§13).

## 4. Имена: переменные, массивы, формулы, функции

Источники имён (`equ_parser.py`, строки 22–61):

| категория | откуда берётся | как индексируется в Fortran |
|---|---|---|
| скаляры (`variables`, `variables_x`, `constants`, `control`, `internInt`, `internDbl`) | `astra_variables.json` (группы описаны в [glossary.md](glossary.md)) | без индекса |
| профили (`profiles`, `profiles_x`, `strahl`) | `astra_variables.json` | `VAR(J)` в радиальном цикле |
| формулы | файлы `fml/<имя>` без `~` и `.` в имени | текст формулы вставляется **перед** оператором (`parse_as.py::fml_fnc`, строки 421–439) |
| функции | файлы `fnc/<имя>.f`, `.f90`, `.F`, `.c`, `.C` (последний символ `f 0 F c C`) | `<ИМЯ>R(RHO(J))`; объявления — в `src/tmp/declar.fnc` |

Скалярные `X`-переменные (`IPLX`, `ZRD1X` …) берутся из эксперимента. Профильные `X`-массивы (`TEX`, `NEX`, `CAR1X` … `CAR64X`, `F0X` … `F9X`, `MUX`, `CUX` …) — тоже. В A7 было только `CAR1X`–`CAR16X` и `F1X`–`F3X` (`section5.tex`, Table 5.4). В A8 есть `CAR1X`–`CAR64X` и `F0X`–`F9X` (`astra_variables.json`, группа `profiles_x`, 106 имён).

Ограничения (`config.py`): `NFML=500` формул и функций (проверяется в `equ_parser.py`, строки 34–37), `NRW=128` радиальных выводов, `NCVA=512`, `NSBMX=20`. Если существуют одновременно функция `fnc/X` и формула `fml/XR`, это фатальная ошибка: `Error! X is both fml and fnc` (`equ_parser.py`, строки 39–42).

Неявная типизация. `write_declar_fml` (`parse_as.py`, строки 241–286) объявляет локальные переменные формул по правилу Fortran: имена на `I`–`N` целые, остальные `double precision`. Это существенно для формул вроде `ZICAR`, `T`, `Z`.

## 5. Выражения и их перевод в Fortran

Трансляцию выполняет `parse_as.py::LINE2FOR` (строки 502–545). Строка режется по `+ - * / ( ) , =` (`rec_split`), скобки разбираются рекурсивно изнутри наружу (`recParse`, `ParseBracket`, `parse_pieces`, строки 323–419). Числа приводятся к двойной точности (`format_number`: `0.5` → `5.D-1`, `1.e-3` → `1.d-3`). Перед разбиением экспоненциальная запись защищается заменой `E-`/`E+` на `$`/`#` (строки 510–525).

Таблица получена прогоном `LINE2FOR` в копии парсера под `/tmp/astra_try`:

| модель | Fortran |
|---|---|
| `CF1=TE(0.3)` | `CF1=RADIAL(TE,RFA(3.D-1))` |
| `CF1=TE(AFX(0.5))` | `CF1=RADIAL(TE,RFA(AFX(5.D-1)))` |
| `CF1=TEB` / `CF1=TEC` | `CF1=TE(NA1)` / `CF1=TE(1)` |
| `CAR1=QETOT` (функция) | `CAR1(J)=QETOTR(RHO(J))` |
| `CF1=QETOTB` | `CF1=QETOTR(ROC)` |
| `CF1=QETOT(0.5)` | `CF1=QETOTR(5.D-1)` |
| `CAR2=VINT(CAR1)` | `CAR2(J)=VINT(CAR1,J*HRO)` |
| `CF1=VINT(CAR1B)` / `IINT(CUBMB)` | `VINT(CAR1,ROC)` / `IINT(CUBM,ROC)` |
| `CAR1=GRAD(TE)` | `enddo` / `do J=1, NA1` / `CAR1(J)=GRAD(TE,J)` |
| `CAR1=1.e-3*NE*TE` | `CAR1(J)=1.d-3*NE(J)*TE(J)` |
| `CAR5=FPR**2` (`FPR` — формула) | `FPR=1.000001-(RHO(J)/ROC)**2` / `CAR5(J)=FPR**2.D0` |
| `CAR4=TE(J+1)` | `CAR4(J)=RADIAL(TE,RFA(J+1.D0))` — индекс **не** поддерживается |
| `CF2=ZRD1X`, `CF5=NA1`, `CF6=TIME` | без изменений |

Правила, видимые в коде (`parse_as.py::indiciseVar`, строки 156–181; `parse_pieces`):

- **Профиль без аргумента** в присваивании профилю даёт `VAR(J)`, оператор попадает в радиальный цикл `do jdetv=1,NA1`.
- **Суффикс `B`** означает значение на границе: `VAR(NA1)` или `fncR(ROC)`. **Суффикс `C`** означает значение в центре: `VAR(1)` или `fncR(0.d0)`. Суффикс распознаётся, если имя без последней буквы — профиль или функция.
- **Профиль с аргументом** `TE(x)` даёт `RADIAL(TE, RFA(x))`, где `x` — малый радиус `a`. Функция с аргументом `F(x)` даёт `FR(x)`, где аргумент передаётся как есть. По коду аргумент функции — это `ρ`, а не `a`. Совпадает ли смысл аргумента у профилей и функций — не проверено. В A8 рекомендуется только форма `TE(AFX(x))` (`docu/addendum_A8.tex`, «Changes in the model file syntax», п. 2).
- **Интегралы** `VINT`, `IINT`, `LININT` от профиля дают `VINT(X, J*HRO)`, от `XB` дают `VINT(X, ROC)`.
- **Функции с аргументами, передаваемые как есть:** `RFVAL RFVEX RFVIN AFVAL AFVEX AFVIN AFX ATX ASTEP RSTEP XSTEP GRAD GRADS RADIAL AFR ATR V_95_POS RFMAX RFMIN FRMAX FRMIN` (`parse_pieces`).
- **`GRAD`/`GRADS`** разрывают радиальный цикл: перед оператором вставляется `enddo` / `do J=1,NA1`, чтобы предыдущие профили были посчитаны на всей сетке.
- **Формулы** раскрываются рекурсивно: `fml.py::IncludedFML` следует строкам `include '…fml/<имя>'`. Текст формулы приводится к верхнему регистру, строки с `!` и `INCLUDE` убираются (`fml.py::textBlock`). В A8 формулу и `GRAD` нельзя писать в одной строке (`addendum_A8.tex`, п. 4).
- **LHS из `internInt`** (`NEQUIL`, `IPEQL`, `NA1` …): строка возвращается без трансляции (`LINE2FOR`, строки 528–529).
- **Кавычки** `"…"` удаляются из результата (строка 540). Так работают вставки `"WORK(J,3)"` в выводе.

**Скалярные присваивания и приоритет данных.** Для скаляров из `variables` генератор `DETVAR` пишет (`code_gen.py`, пример — в `~/astra/a8/src/tmp/detvar.f90`):

```fortran
IFDFVX(19) = max(IFDFVX(19), 2)
if (IFDFVX(19) <= 2) IPL=IPLX
```

`IFDFVX` — источник значения: 0 — эксперимент без времени, 1 — эксперимент с временем, 2 — модель, 3 — клавиатура, 4 — запрещено менять (`src/astra/set_x_data.f90::set_x_scalars`). Значение из модели перекрывает эксперимент, но не перекрывает ввод с клавиатуры. Это аналог иерархии A7 (`section5.tex`, Table 5.1).

## 6. Уравнения переноса

### 6.1. Тип уравнения

`equ_parser.py`, строки 175–181, задаёт `LEQ`, а `eqns_init.py` и `eqns.py` генерируют код. Таблица проверена прогоном парсера под `/tmp/astra_try` (модель из одной строки `TE:<тип>`):

| запись | `LEQ(2)` в `ininam.f90` | код в `eqns_inc.f90` |
|---|---|---|
| `TE:EQ` / `TE:Equation` / `TE:EQ[…]` | 1 | «Electron temperature equation» (`teeqn`) |
| `TE:AS` | 0 | «Electron temperature assignment» |
| `TE:;` (пустой тип) | 1 | **нет ни уравнения, ни присваивания**, профиль заморожен на начальном значении |
| `TE:E`, `TE:A` (сокращения A7) | −1 | **нет кода**, без предупреждения |
| `TE:FU` | 2 | по коду ветки нет, результат не проверен |
| строки `TE:` нет | −1 | присваивание `TE=TEX` при старте (`Missing`) |

Причина в разных проверках. `LEQ` считается по `val[:2] in ('', 'EQ')`, `val == 'AS'`, `val == 'FU'`. Генераторы проверяют `assign_type[:2] == 'EQ'` и `== 'AS'` (`eqns.py::pre_eqn`, строки 25–36; `eqns_init.py`). Для тока `CU:;` даёт `LEQ(4)=1`, но в коде только «CU adjustment» (`cuasn`), а не уравнение (проверено). В A7 `CU:;` означало решение уравнения (`section5.tex`). **В A8 тип нужно писать полностью: `EQ` или `AS`.**

Неявная связка: `TE*:EQ` / `TI*:EQ` даёт `assign_d = 'implicit_EQ'`, и `eqns_init` вызывает `tetieqn` вместо `teeqn`+`tieqn` (`eqns_init.py`; `eqns.py`, строки 645–844). Этой формы в A7 не было.

### 6.2. Что генерируется для уравнения

`pre_eqn` (`eqns.py`, строки 17–51) открывает цикл `do J=1, NA1`. В цикле вычисляются все коэффициенты уравнения, неопределённые источники обнуляются (`PE(J) = 0.`), и копируется полный поток `PETOT(J)=PE(J)`. Затем `teeqn`/`tieqn`/`neeqn` (строки 251–456) собирают коэффициенты прогонки `YWA`…`YWD` и вызывают решатель (см. [05-solver-core.md](05-solver-core.md)). Для `AS` генерируется присваивание профиля из правой части и пересчёт потоков (`QE`, `QI`, `QN`).

Ток (`cueqn`, строки 456–513): если задан `IPL`, полный ток задаётся как граничное условие; варианты с `UEXT`/`LEXT` описывают внешнюю цепь. `CU:AS` или отсутствие `EQ` вызывает `cuasn` (строки 182–224): ток задаётся профилем `CU` или `MU`. **Если в модели нет ни `CU`, ни `MU`, парсер падает** с `UnboundLocalError: … 'cu_as'` (`eqns_init.py`, строка 48; проверено).

Примеси и прочие уравнения `F0`…`F9` (`fjeqn`, строки 515–588) вызывают `RUNEQ`. Вращение `UPAR` обрабатывает `upeqn` (строки 590–643).

### 6.3. Граничные условия и `EQ[...]`

`rhoBC` (`eqns.py`, строки 54–81) разбирает квадратные скобки:

| запись | граница уравнения `ND1` | окно по времени |
|---|---|---|
| `TE:EQ` | `ND1 = NA1` | — |
| `TE:EQ[x]` | `NODE(RFA(x))`, `x` — `a` в метрах | всегда |
| `TE:EQ[1,x]` | `NODE(RFAN(x))`, `x` — нормированный `ρ` | всегда |
| `TE:EQ[2,x]` | `NODE(x*ROC)` | всегда |
| `TE:EQ[1,x,t1]` | как выше | с `t1` |
| `TE:EQ[1,x,t1,t2]` | как выше | `[t1,t2]` |

Окно по времени — расширение A8 (`bnd_init`, строки 120–141). При `TIME<t1` ставится `ND1=1`, профиль целиком задаётся внешним выражением. При `TIME>t2` ставится `ND1=1`, и профиль замораживается на `VARO`. Первый аргумент разбирается как `int(tmp[0])`, то есть берётся только первый символ. Значения, отличные от 1 и 2, вызывают `sys.exit()`.

Тип граничного условия (`bnd_init`, строки 146–177):

| что задано в модели | результат |
|---|---|
| ничего (`TEB`, `QEB`, `QETB` нет) | `bctype=1`. Если и внешнего выражения `TE=…` нет, выдаётся предупреждение и ставится `TE(ND1:NA1)=TEX(ND1:NA1)` |
| `TEB=…` | `bctype=1`, `TE(ND1)=…`. Если внешнего выражения нет, между `ND1` и сепаратрисой ставится линейная интерполяция (`linInterp2sep`) |
| `QEB=…` или `QETB=…` | `bctype=2`, задан поток `QE(ND1)` (для `QETB` поток равен `TE(ND1)*(…)`) |
| больше одного | `ValueError('Too many boundary conditions …')` |

Правило A8 (`addendum_A8.tex`, п. 3): начальное условие `NE=NEX` означает также `NEB=NEX` на границе `EQ[...]`, если `NEB` не задан явно. В коде это видно в цикле `do j=ND1,NA1` с внешним выражением.

## 7. Подпрограммы

### 7.1. Синтаксис

```text
NAME(arg1, arg2, ...)<:dt:tmin:tmax:key
```

Разбор выполняет `parse_as.py::parse_sbr` (строки 441–498). Имя берётся до `(` и `:`, символы `<` и `>` по краям срезаются. Аргументы передаются в `call NAME(...)` без трансляции, но в верхнем регистре.

Число двоеточий определяет, какие параметры заданы:

| двоеточий | поля | значения по умолчанию, которые ставит парсер |
|---|---|---|
| 4 | `dt:tmin:tmax:key` | пустое поле не записывается, остаётся умолчание ядра |
| 3 | `dt:tmin:tmax` | `key=''` |
| 2 | `dt:tmin` | `tmax=1000.` |
| 1 | `dt` | `tmin=0.`, `tmax=1000.` |

Умолчания ядра (`src/astra/scalars.f90`, строки 157–160): `DTEQ(1,:)=0`, `DTEQ(2,:)=-99999`, `DTEQ(3,:)=99999`, `DTEQ(4,:)=-1`. Отсюда разница в поведении: **`NEUT:.05` не вызывается при `t<0` и `t>1000`, а `NEUT:.05:::` вызывается всегда** (по коду, на счёте не проверено).

Буква `key` переводится в номер `ord(key)-ord('a')+1` (строка 496). Подпрограмма вызывается при нажатии этой клавиши в GUI: `ABS(KEY - DTEQ(4,j)) < 0.1`.

Если поле — число, оно идёт в `ininam.f90` (`DTEQ(1,1) = 5.d-2`). Если поле — выражение, оно идёт в `detvar.f90` и пересчитывается каждый шаг (`code_gen.py`, строки 60–69).

Условие вызова (`parse_as.py::sbr_header`, строки 288–297; сгенерировано в `eqns_inc.f90`):

```fortran
IFSUB = 0
if (KEY /= 0 .and. ABS(KEY - DTEQ(4,1)) < 0.1) IFSUB = 1
if (TIME >= DTEQ(2, 1) .and. TIME <= DTEQ(3, 1) .and. TIME - TEQ(1) + 1.E-7 - DTEQ(1, 1) > 0.0) IFSUB = 1
```

Каждый вызов обёрнут замером `wallTime_sbr`/`cpuTime_sbr` (`write_sbr`, строки 299–321). Для `STRAHL` есть особая ветка.

### 7.2. Место вызова (`locsbr`)

| метка | `locsbr` | файл | когда вызывается |
|---|---|---|---|
| нет | 0 | `init_converge_step.f90` и `eqns_inc.f90` | при старте (цикл сходимости) и на каждом шаге в блоке уравнений |
| `<` в строке | −1 | `detvar.f90` | в `DETVAR` при `IPART==2`, до равновесия |
| `<<` | задумано −2 | — | **недостижимо**: после `if '<<'` идёт `if '<'` без `elif`, поэтому всегда −1 (строки 473–476; проверено на `NEUT<<:0.1:::` → `detvar.f90`) |
| `>` в строке, или имя `MIXINT`/`MIXEXT`/`TSCTRL` | 1 | `postep.f90` | после шага (`stepup.F90`, `POSTEP`, строка ~240) |

Проверки `'<' in line` и `'>' in line` ищут символ по всей строке, включая аргументы. Например, `=>` в аргументе отправит вызов в `postep` (по коду).

### 7.3. Именованные аргументы

`RABBIT(pRF_MW=1.d0, …)<:.05:::;` попадает одновременно в ветку `=` (как «присваивание» `RABBIT(PRF_MW = …`, которое `code_gen` игнорирует, потому что LHS не переменная) и в ветку `:` (подпрограмма). Аргументы `PRF_MW=1.D0` передаются в Fortran как ключевые, то есть у подпрограммы должен быть явный интерфейс (модуль). Так сделано в `a2rabbit` (`detvar.f90`: `use a2rabbit, only: rabbit`).

## 8. Вывод

Подробно графика описана в [08-graphics-and-output.md](08-graphics-and-output.md). Здесь — только синтаксис.

**Радиальный вывод** (`equ_parser.py`, строки 235–262):

```text
Name\Expr                 ! профиль
Name\Expr\Scale           ! с масштабом; отрицательный масштаб группирует окна
Name\Expr\\XDATA          ! с экспериментальными точками
Name\Expr\\XDATA\Scale
\                         ! пустое окно
```

Из `\\…` берётся только **последняя** ссылка: `sgr.split('\\\\')[-1]`. В A7 можно было писать несколько. Число окон ограничено 128 (`NRW`), лишние отбрасываются с предупреждением. Выражение переводится через `LINE2FOR` в `astra_out.f90` (`ROUT(…)`).

**Временной вывод** (строки 220–231):

```text
Name_Expr
Name_Expr_Scale
```

Выражение режется по `_`, поэтому `_` в самом выражении недопустим. Строка с `_` без `:`, `\`, `=` всегда считается выводом, даже если это задумывалось как комментарий.

## 9. Отличия A8 от A7

Сводка по `docu/addendum_A8.tex` (февраль 2024, «Changes in the model file syntax»), `docu/section5.tex` (руководство A7) и коду:

| тема | A7 | A8 |
|---|---|---|
| транслятор | собственный построитель (в репозитории нет) | `pyparse/*.py` → `src/tmp/*.f90` |
| функции `fnc` | как функция `ρ`: `QETOT(r)` | только как массив или граничное значение `QETOTB` (addendum, п. 1); по коду `F(x)` → `FR(x)` всё ещё генерируется |
| значение профиля в точке | `TE(0.3)` | только `TE(AFX(0.3))` (п. 2) |
| граница уравнения при `NE=NEX` | — | подразумевается `NEB=NEX` на `EQ[...]` (п. 3) |
| формула и `GRAD` | в одной строке | нельзя (п. 4) |
| тип уравнения | `TE:;`, `TE:E`, `TE:A`, полные слова | только `EQ…` и `AS` (§6.1) |
| окно по времени в `EQ[...]` | нет | `EQ[1,x,t1,t2]` |
| неявная связка TE–TI | нет | `TE*:EQ` |
| `X`-массивы | `CAR1X`–`CAR16X`, `F1X`–`F3X` | `CAR1X`–`CAR64X`, `F0X`–`F9X` |
| `\\XDATA` | несколько | одна (последняя) |
| `<<` (до равновесия) | — | недостижимо (§7.2) |
| препроцессор `-i template` | — | описан в `addendum_A8.tex`, но опции `-i` в `repo_user_area/exe/as_exe.py` нет |
| подпрограммы `STARR`, `STVAR`, `TSCTRL`, `MINMAX`, `SETNAV`, `GNXSRC` | были | в `sbr/` нет ([09-user-area-catalog.md](09-user-area-catalog.md)) |

Перенос моделей A7 на A8 по пунктам:

1. Заменить сокращённые и пустые типы на `EQ` / `AS`.
2. Заменить `F(r)` и `TE(r)` на `FB` или `TE(AFX(…))`.
3. Разнести формулы и `GRAD` по разным строкам.
4. Оставить одну `\\XDATA` на окно.
5. Задать ток (`CU` или `MU`), иначе парсер падает.
6. Проверить, что все вызываемые подпрограммы есть в `sbr/`.

## 10. Модули `pyparse` по отдельности

| модуль | строк | назначение |
|---|---|---|
| `parser_main.py` | 58 | точка входа: `astra_parser(f_equ, f_exp)` → `EQU_PARSER`, `EXP_PARSER`, предупреждения, `nr_x_max`, `CODE_GEN`; `write_tmp` пишет в `./src/tmp` |
| `config.py` | 63 | пути (`awd` — родитель `pyparse`, `fml_dir`, `fnc_dir`), лимиты, `eqn_list`, `coeff_d`, `flux_d`, `bnd_d`, подписи |
| `parse_as.py` | 545 | предобработка, трансляция выражений `LINE2FOR`, подпрограммы `parse_sbr`/`sbr_header`/`write_sbr`, объявления, запись Fortran |
| `equ_parser.py` | 260 | `EQU_PARSER`: списки имён, классификация строк, `LEQ`, начальные выражения `init_d`, вывод |
| `code_gen.py` | 371 | `CODE_GEN`: собирает тексты всех `src/tmp/*` |
| `eqns_init.py` | 89 | порядок блоков уравнений в `eqns_inc` и `init_converge_step` |
| `eqns.py` | 844 | генераторы уравнений и граничных условий: `pre_eqn`, `rhoBC`, `bnd_init`, `neeqn`, `teeqn`, `tieqn`, `tetieqn`, `cueqn`, `cuasn`, `fjeqn`, `upeqn` |
| `const_text.py` | 971 | шаблоны Fortran: заголовки `DETVAR`, `INIVAR`, `ININAM`, `POSTEP`, `SET_GRAPH_NAMES`, блоки уравнений, `RADOUT`/`TIMOUT` |
| `fml.py` | 79 | рекурсивное раскрытие формул (`IncludedFML`, `FML.textBlock`) |
| `rw_for.py` | 30 | ввод-вывод массивов Fortran-формата (`'%13.6E'`, 6 в строке), разбор слипшихся чисел `1.0-2.0` |
| `test_parse.py` | 25 | `python3 test_parse.py -s "HE = …"` печатает результат `LINE2FOR` |
| `exp_parser.py` | 305 | разбор файла эксперимента на этапе сборки (см. [04-experiment-data.md](04-experiment-data.md)) |
| `ufiles.py` | 345 | чтение и запись U-файлов (Tardini, v0.3.0, 2022) |
| `repo_user_area/pyparse/iondens_ass.py` | — | `NIAS.iondensassign`: `NI=F1+…+F9`, включается строкой с `@!@` |

**`code_gen.CODE_GEN`** собирает тексты по частям:

- **`ininam`**: `AWD`, `LEQ(1..16)`, `sbr_name(j)`, числовые `DTEQ`.
- **`associate_pointers`**: указатели на массивы `constValues`, `varValues`, `varxValues`, `controlValues`, `internInt`, `internDbl`, `profiles_x`, `profiles`, а также `n_sbr` и `nr_x_max`.
- **`detvar`**: скаляры; подпрограммы с `<`; радиальный цикл `do jdetv=1,NA1; J=jdetv`; вставки в кавычках; хвост для `IPEQL==3`.
- **`inivar`**: `F0`…`F9` по умолчанию 1; начальные `NE`/`TE`/`TI`/`UPAR` (из `init_d` или `X`-массивов); ограничение `NI`; `CC=CCSP` по умолчанию; `MU`/`CU`.
- **`postep`**: подпрограммы с `>`.
- **`set_graph_names`**, **`astra_out`**: вывод.
- **`model.txt`**: читаемая сводка.
- **`eqns_inc`**, **`init_converge_step`**: через `eqns_init`.

Ошибка трансляции. `LINE2FOR` ловит исключение `recParse` и печатает `Error in EQU-file line:` с traceback, но после этого обращается к неопределённой `line_out`. По коду парсер затем падает с `UnboundLocalError`, а не продолжает работу.

## 11. Сгенерированные файлы и их место в счёте

`parser_main.write_tmp` пишет в `src/tmp/`. `Makefile` компилирует `associate_pointers`, `inivar`, `detvar`, `postep`, `ininam`, `init_converge_step`, `eqns_inc`, а при сборке с X11 ещё `astra_out` и `set_graph_names`. Флаг `-I$(AWD)` позволяет `include 'src/tmp/declar.fml'`.

| файл | содержимое | кто вызывает (`src/astra/astra_main.F90`, `stepup.F90`) |
|---|---|---|
| `ininam.f90` | `LEQ`, имена и интервалы подпрограмм | `ininam` при старте, до `readInput` |
| `associate_pointers.f90` | связь имён `astra_variables.json` с памятью | первым после `read_metadata` |
| `inivar.f90` | начальные профили | `INIVAR` при старте и в цикле сходимости |
| `detvar.f90` | скаляры и профили модели, подпрограммы `<` | `DETVAR` при старте, в цикле сходимости и на каждом шаге |
| `init_converge_step.f90` | подпрограммы и уравнения для итераций до `ATREQ` | `INIT_CONVERGE_STEP` в стартовом цикле |
| `eqns_inc.f90` | подпрограммы без метки, уравнения | шаг по времени |
| `postep.f90` | подпрограммы `>` | `POSTEP` в `stepup.F90` |
| `set_graph_names.f90`, `astra_out.f90` | имена и выражения вывода | графика, `write_ajson` |
| `declar.fml`, `declar.fnc` | объявления локальных переменных формул и функций `…R` | `include` в генерируемых файлах |
| `model.txt` | сводка выводов | человек |

Стартовая последовательность (`astra_main.F90`): `read_metadata` → `associate_pointers` → `ininam` → `readInput` → `astra_assignments` → `set_x_arrays(1)` → `INIVAR` → `SETVAR` → `DETVAR` → … затем цикл до `IFTREQ(ATREQ)`: `set_x_scalars`, `DETVAR`, `DEFARR`, `set_x_arrays(1)`, `INIVAR`, `INIT_CONVERGE_STEP`, `METRIC`. На шаге времени (`stepup.F90`) идут `set_x_scalars` (строка ~65), `set_x_arrays(2)` (~119), уравнения и `POSTEP` (~240). Детали — в [05-solver-core.md](05-solver-core.md).

## 12. Разбор примера: `flux_feqis` → `src/tmp`

Модель `repo_user_area/equ/flux_feqis` (109 строк) запускается с `aug34954` (`tmp/astra.nml`: 4.0–5.0 с, `BGD`, `machine='aug'`). В `~/astra/a8/equ/flux_feqis` установщик закомментировал строки `TORBA` и `RABBIT` (53–54). Это единственное отличие от исходника.

| строка модели | результат в `~/astra/a8/src/tmp` |
|---|---|
| `NEQUIL=91; MEQUIL=90; IPEQL=5;` | `detvar.f90`: `NEQUIL = 91` … без трансляции (LHS из `internInt`) |
| `atreq = .5;` | `ATREQ=5.D-1` |
| `IPL=IPLX;` | `IFDFVX(19) = max(IFDFVX(19), 2)` / `if (IFDFVX(19) <= 2) IPL=IPLX` |
| `AIM1=12.;` | то же с `IFDFVX(3)` |
| `TAUMIN=0.05;` | `TAUMIN=5.D-2` (переменная `control`, без `IFDFVX`) |
| `ZIM1=ZICAR;` | в радиальном цикле перед оператором вставлен текст формулы `fml/zicar` (`T=TE(J)`, `Z=LOG10(T)`, `IF …`) |
| `NIZ1=CIMP1*NE;` | `NIZ1(J)=CIMP1*NE(J)` |
| `CF15=BETANB;` | `CF15=BETANR(ROC)` |
| `CF6=0.33*VINT(PRADXB)+…` | `CF6=3.3D-1*VINT(PRADX,ROC)+VINT(CAR15,ROC)+VINT(CAR16,ROC)` |
| `CF8=NEC/NEAVB;` | `CF8=NE(1)/NEAVR(ROC)` |
| `CF9=TIX(AFX(0.));` | `CF9=RADIAL(TIX,RFA(AFX(0.D0)))` |
| `NEUT:.05:::;` | `ininam.f90`: `sbr_name(1) = "NEUT"`, `DTEQ(1,1) = 5.d-2`; вызов в `eqns_inc.f90` и `init_converge_step.f90` |
| `TE:AS; TE=TEX;` (то же `TI`, `NE`) | `LEQ(1..3)=0`; в `eqns_inc.f90` блоки «… assignment»; в `inivar.f90` `TE(J)=TEX(J)` |
| `CU:EQ; CC=CNSA; CD=CUBM+CUECR;` | `LEQ(4)=1`, блок уравнения тока; `CC` — раскрытая формула `CNSA` |
| `MU=1./(CAR2X + 0.0001);` | `inivar.f90`: начальный `MU` |
| строка `Equil` отсутствует | `LEQ(5) = -3` (всегда, см. §13) |
| `Profile output` / `Time output` | комментарии (нет спецсимволов) |
| `Tex\TE\\TEX;` … `q\1./MU\\CAR2X;` | `set_graph_names.f90`: `NROUT=32`, `NXOUT=4` (`TEX NEX TIX CAR2X`, окна 1, 5, 9, 14) |
| `QE_QETOTB; …` (18 штук) | `NTOUT=18` |

Размеры сгенерированных файлов в `~/astra/a8/src/tmp` (`wc -l`): `associate_pointers` 1119, `init_converge_step` 315, `eqns_inc` 445, `detvar` 161, `astra_out` 168, `inivar` 92, `set_graph_names` 87, `ininam` 33, `postep` 31 (пустое тело).

Фрагмент `eqns_inc.f90` для `NEUT:.05:::` и присваивания `NE`:

```fortran
IFSUB = 0
if (KEY /= 0 .and. ABS(KEY - DTEQ(4,1)) < 0.1) IFSUB = 1
if (TIME >= DTEQ(2, 1) .and. TIME <= DTEQ(3, 1) .and. TIME - TEQ(1) + 1.E-7 - DTEQ(1, 1) > 0.0) IFSUB = 1
if (IFSUB == 1) then
    TEQ(1) = TIME
    call markloc("subroutine NEUT")
    ...
    call NEUT()
    ...
endif
! **** Density assignment
call markloc("NE assignment")
do J=1, NA1
    SN(J)=SNEBM(J)
    SNN(J) = 0.
    SNTOT(J) = SN(J)
enddo
```

Поля `tmin`/`tmax` в `NEUT:.05:::` пустые, поэтому `DTEQ(2,1)` и `DTEQ(3,1)` остаются −99999/99999 из `scalars.f90`.

## 13. Ловушки и ошибки парсера

| № | ловушка | источник | статус |
|---|---|---|---|
| 1 | `TE:;`, `TE:E`, `TE:A` молча не генерируют уравнение; `CU:;` даёт присваивание тока вместо уравнения | `equ_parser.py`, строки 175–181; `eqns.py::pre_eqn`; `eqns_init.py` | проверено в `/tmp` |
| 2 | Модель без `CU` и `MU` роняет парсер (`UnboundLocalError: cu_as`) | `eqns_init.py`, строка 48 | проверено |
| 3 | `Equil:…` не распознаётся (ключ `Equil` в смешанном регистре против `EQUIL`) и превращается в `call EQUIL()`; `LEQ(5)` всегда −3; таблица `equ_solver` (`EQ0`, `ESC`, `SPIDER` …) недостижима | `equ_parser.py`, строки 144, 164–185; `config.py` | проверено (`call EQUIL()` в `eqns_inc.f90`); исход линковки не проверен |
| 4 | `<<` недостижим, всегда даёт `<` | `parse_as.py`, строки 473–476 | проверено |
| 5 | `NAME:dt` (1 двоеточие) ограничивает вызов окном `[0, 1000]` с, а `NAME:dt:::` не ограничивает | `parse_as.py`, строки 468–472; `scalars.f90`, строки 157–160 | по коду |
| 6 | Аргументы подпрограмм в верхнем регистре, в том числе строковые литералы | `parse_sbr`, строка 448 | по коду |
| 7 | `!` и `;` внутри кавычек режут строку | `equ_prepare`, строки 144–145 | по коду |
| 8 | Явный индекс `TE(J+1)` транслируется в `RADIAL(TE,RFA(J+1.D0))` | `parse_pieces` | проверено |
| 9 | В `model.txt` `var.replace('B','(a)')` заменяет все `B`: `betn_BETANB` → `(a)ETAN(a)` | `code_gen.py`, строка 216 | косметика, видно в `~/astra/a8/src/tmp/model.txt` |
| 10 | При отсутствии подпрограмм комментарий `! **** No external subroutines` пишется без перевода строки, и `call markloc("eqns")` попадает в комментарий | `eqns_init.py`, строка 31 | проверено (`/tmp`), безвредно |
| 11 | Ошибка трансляции выражения → traceback и `UnboundLocalError` | `LINE2FOR`, строки 531–537 | по коду |
| 12 | Строки `_` без `:`, `\`, `=` считаются временным выводом, даже если это задумывалось как текст | `equ_parser.py`, строки 94–102 | по коду |
| 13 | Учёт `X`-массивов (`arname`) берёт из строки не больше одного нового имени (`break`). Это влияет только на предупреждения `parser_main` | `equ_parser.py`, строки 211–216 | по коду |
| 14 | Файл эксперимента не входит в зависимости `make`, поэтому `nr_x_max` может устареть | `repo_user_area/exe/Makefile` | по коду; последствия не проверены |
| 15 | Предупреждение `nrho>500` в эксперименте; лимит `NRD` проверяется Fortran (`astra_assignments`: `NA1>NRD` фатально) | `parser_main.py`; `set_x_data.f90` | по коду |

**Патч U-строк A7.** Установщик заменяет `parse_u_line` в `~/astra/a8/pyparse/exp_parser.py` на регулярное выражение. Так парсер сборки принимает строки `TEX<tab>U-file:…<tab>factor:…` и пропускает строки-комментарии с `U-file` ([DETAILED.md §3.8](../../installer/docs/DETAILED.md)). Патч касается **только** этапа сборки: во время счёта Fortran молча игнорирует такие строки с табуляцией (проверено, см. [04-experiment-data.md](04-experiment-data.md), §9). В исходниках `astra-src` патча нет.

## 14. Как проверять модель без сборки

Генерацию можно запустить в отдельной копии, не трогая рабочее дерево:

```bash
mkdir -p /tmp/astra_try && cd /tmp/astra_try
cp -r ~/astra/a8/pyparse ~/astra/a8/fml ~/astra/a8/fnc ~/astra/a8/astra_variables.json .
mkdir -p equ exp src/tmp
cp ~/astra/a8/equ/flux_feqis equ/; cp ~/astra/a8/exp/aug34954 exp/
python3 pyparse/parser_main.py -equ equ/flux_feqis -exp exp/aug34954
ls src/tmp; grep LEQ src/tmp/ininam.f90
```

`config.awd` — родитель каталога `pyparse`, поэтому копия ищет `fml/`, `fnc/` и `astra_variables.json` рядом с собой. `write_tmp` пишет в `./src/tmp` относительно текущего каталога.

Отдельное выражение проверяется так:

```bash
cd /tmp/astra_try/pyparse && python3 test_parse.py -s "CAR1=VINT(PE)"
```

Минимальная проверка новой модели: в `ininam.f90` нужные `LEQ(..)=1` или `0`; в `eqns_inc.f90` есть заголовки «… equation» для каждого уравнения; `sbr_name(..)` совпадают с файлами в `sbr/`; в выводе парсера нет строк `Error` и `Missing`.

## 15. Где смотреть в исходниках

| файл | что там |
|---|---|
| `pyparse/parser_main.py` | точка входа, `nr_x_max`, предупреждения, запись `src/tmp` |
| `pyparse/config.py` | `eqn_list`, `coeff_d`, `flux_d`, `bnd_d`, лимиты |
| `pyparse/parse_as.py` | `equ_prepare` (118), `indiciseVar` (156), `sbr_header` (288), `write_sbr` (299), `parse_pieces` (354), `fml_fnc` (421), `parse_sbr` (441), `LINE2FOR` (502) |
| `pyparse/equ_parser.py` | классификация строк (90–152), `LEQ` (157–185), вывод (220–262) |
| `pyparse/eqns.py` | `pre_eqn` (17), `rhoBC` (54), `bnd_init` (96), `cuasn` (182), `neeqn` (251), `tieqn` (300), `teeqn` (377), `cueqn` (456), `fjeqn` (515), `upeqn` (590), `tetieqn` (645) |
| `pyparse/eqns_init.py` | порядок блоков уравнений |
| `pyparse/code_gen.py` | генерация всех файлов `src/tmp` |
| `pyparse/const_text.py` | шаблоны Fortran |
| `pyparse/fml.py`, `pyparse/rw_for.py`, `pyparse/test_parse.py` | формулы, ввод-вывод массивов, отладка выражений |
| `repo_user_area/pyparse/iondens_ass.py` | `NI=F1+…+F9` по `@!@` |
| `astra_variables.json` | все имена переменных по группам |
| `repo_user_area/exe/Makefile` | `PY_TMP_GEN`, `TMP_SRC`, зависимости |
| `repo_user_area/exe/as_exe.py` | опции командной строки, `tmp/astra.nml` |
| `src/astra/astra_main.F90`, `src/astra/stepup.F90` | порядок вызова сгенерированных подпрограмм |
| `src/astra/scalars.f90` (157–160) | умолчания `DTEQ` |
| `src/astra/set_x_data.f90` | `IFDFVX`, `set_x_scalars`, `astra_assignments` |
| `docu/addendum_A8.tex` | изменения синтаксиса A8 |
| `docu/section5.tex`, `docu/section4.tex` | руководство A7: синтаксис, иерархия данных, встроенные функции |
| `repo_user_area/equ/flux_feqis`; `~/astra/a8/src/tmp/*` | разобранный пример |
