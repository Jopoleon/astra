# 08. Графика, вывод результатов и регрессии

**О чём глава.** Как ASTRA 8 показывает расчёт и что от него остаётся на диске. Разобраны: X11-интерфейс (`src/graph`, `lib/libgraph.a`), его кнопки, клавиши и режимы; объявления экранного вывода в модели (`Profile output`, `Time output`) и почему они не попадают в CDF; JSON-вывод ядра (`write_ajson`); склейка JSON в NetCDF (`postProcess/json2cdf.py`); структура и именование `ncdf_out/<exp><equ>.CDF` с проверенным на реальном файле примером чтения на Python; прочие скрипты `postProcess/`; эталоны `regressions/`, `compReg.sh`, `compareRegressions.py`, допуск сравнения и известные расхождения gfortran с Intel. Состояние кода — `astra-src` `67281112` (2026-10-01); все утверждения о формате CDF проверены на `~/astra/a8/ncdf_out/aug34954flux_feqis.CDF` (расчёт `flux_feqis`/`aug34954`, `-s 4 -e 5`).

Что **не** повторяется здесь:

- запуск, опции `as_exe`, GUI/batch на уровне конвейера — [02-build-and-run.md](02-build-and-run.md);
- синтаксис разделов модели — [03-model-language.md](03-model-language.md);
- экспериментальные данные и `...X`-массивы — [04-experiment-data.md](04-experiment-data.md);
- величины равновесия в группах `equil_*` — [06-equilibrium.md](06-equilibrium.md);
- откуда берутся `PEBM`, `PEECR` и т. п. в регрессиях — [07-heating-and-external-modules.md](07-heating-and-external-modules.md);
- термины — [glossary.md](glossary.md).

## Содержание

1. [Каналы вывода одним взглядом](#1-каналы-вывода-одним-взглядом)
2. [X11-интерфейс: устройство](#2-x11-интерфейс-устройство)
3. [Кнопки и клавиши](#3-кнопки-и-клавиши)
4. [Режимы графики](#4-режимы-графики)
5. [Экранный вывод в модели](#5-экранный-вывод-в-модели)
6. [JSON-вывод ядра](#6-json-вывод-ядра)
7. [`json2cdf.py`: склейка в NetCDF](#7-json2cdfpy-склейка-в-netcdf)
8. [Формат CDF](#8-формат-cdf)
9. [Чтение CDF на Python](#9-чтение-cdf-на-python)
10. [Рестарт из JSON](#10-рестарт-из-json)
11. [Скрипты `postProcess/` и `plot_compare.py`](#11-скрипты-postprocess-и-plot_comparepy)
12. [Регрессии: `regressions/`, `compReg.sh`, `compareRegressions.py`](#12-регрессии-regressions-compregsh-compareregressionspy)
13. [Документация A7 против реальности A8](#13-документация-a7-против-реальности-a8)
14. [Известные дефекты и подводные камни](#14-известные-дефекты-и-подводные-камни)
15. [Где смотреть в исходниках](#15-где-смотреть-в-исходниках)

---

## 1. Каналы вывода одним взглядом

| Канал | Кто пишет | Куда | Когда | Содержимое |
|---|---|---|---|---|
| JSON-снимки | `write_ajson` (`src/astra/json_rw.f90`) | `ncdf_out/<exp><equ>-<N>.json` | каждые `DPOUT` секунд модельного времени | все внутренние переменные ядра, профили, равновесие |
| NetCDF | `postProcess/json2cdf.py` (вызывается из `exe/as_exe.py`) | `ncdf_out/<exp><equ>.CDF` | после завершения расчёта | склейка всех JSON по времени |
| X11-окно | `src/graph/*` | экран | в GUI-сборке (без `-batch`) | экранный вывод модели, профили, равновесие |
| Текстовый дамп экрана | клавиша `F` (`writeData`) | `dat/<exp>-<equ>-<k>.dat` | по запросу | текущие значения **экранных** выражений модели |
| PostScript | клавиши `Q`/`G` | `dat/<exp>-<equ>-<k>.ps` | по запросу | текущий кадр |
| Лог параметров | клавиша `I` | `equ/log/<equ>` | при входе и по запросу | переменные, константы, управляющие параметры |
| EQDSK | `EQDSK(coco)` в модели (`sbr/a2eqdsk.f90`) | `ncdf_out/<exp><equ><время>.eqdsk` | по расписанию вызова | равновесие в формате G-EQDSK |
| Модульные файлы | NBI, TORIC, STRAHL и др. | `dat/`, `toric/`, `strahl/<exp>_<equ>/` | по вызову модулей | см. [07](07-heating-and-external-modules.md) |
| Терминал | всё ядро | stdout | всегда | диагностика; успешный конец — `>>> ASTRA normal exit >>>` |

`k` в именах `.dat`/`.ps` — первый свободный номер 1, 2, 3, … (функция `set_filename`, `src/graph/gui_interaction.f90` ≈1334–1356).

Главное различие: **CDF содержит только внутренние переменные ядра**, а не то, что модель выводит на экран. Подробно — §5.

---

## 2. X11-интерфейс: устройство

### 2.1. Сборка

GUI существует только при сборке с X11. В `repo_user_area/exe/Makefile` (≈52–72): `ifdef X11` → поиск `X11_LIB` (проба multiarch-каталогов), флаг `-DX11`, линковка `-lgraph`, и к сгенерированным исходникам добавляются:

- `src/tmp/astra_out.f90` — подпрограммы `RADOUT`/`TIMOUT`, вычисляющие экранные выражения модели (`ROUT(J,k)=<выражение>`);
- `src/tmp/set_graph_names.f90` — имена и масштабы (`NAMER`, `SCALER`, `NAMET`, `NTOUT`, `NROUT`, `NXOUT`; для `flux_feqis`: `NTOUT=18`, `NROUT=32`, `NXOUT=4`).

Генерирует их `pyparse/code_gen.py` (≈311–358) по шаблонам `pyparse/const_text.py` (`class RADOUT`, `class TIMOUT`, ≈586–624). Следствие: **в batch-сборке экранные выражения модели вообще не компилируются**.

### 2.2. Компоненты `src/graph`

| Файл | Роль |
|---|---|
| `Astra2XW.c` | низкоуровневая X11-графика: `initvm_`, `drawvm_`, `textvm_`, `resizewindow_`, `close_screen_`, вывод PostScript `psopen_` |
| `dialogx.c` | меню из 24 кнопок (`taskmenu`), диалоговые таблицы `MENUTABLE` с заголовками `Variable control`, `Constant control`, `Times and grids`, `Sequence control`, `Time interval`; заголовки групп меню `Graphic mode`, `Select`, `Control`, `In`, `Out`, `Status` |
| `graph_utils.f90` | `gui_init` (окно 660×550, умножается на `-resize`, заголовок `Per aspera ad ASTRA`, версия из `exe/version`), `plotMode`, `set_plot_area`, `show_plots`, `plot_lcfs`, `plot_wall`, `plot_flux_surfaces`; лимиты `nplots_max=32`, `NRW=128`, `NTIMES=1024`, `plot_modes=9` (строка 10) |
| `gui_interaction.f90` | обработчик клавиш `if_key`, `writeData`, сохранение лога, PostScript |

### 2.3. Жизненный цикл

`src/astra/astra_main.F90`: при `TASK≠BGD` (и `#ifdef X11`) вызывается `set_graph_names` и `gui_init`. Во время начальных итераций и каждого шага вызывается `if_key` — окно перерисовывается и опрашиваются клавиши. В GUI-режиме расчёт **стартует сам** и при `TIME ≥ TEND` завершается вместе с окном (`STOP` после `CPU_report`). Просмотрщика результатов после конца расчёта в A8 нет (§13). Как выглядит окно тестового прогона — `installer/docs/DETAILED.md`, раздел 5.2.

Аргумент `-tpause T` у `as_exe` ставит автоматическую паузу при `TIME=T` (по умолчанию `1e4`, т. е. никогда).

---

## 3. Кнопки и клавиши

### 3.1. Кнопки меню

Из `src/graph/dialogx.c` (≈140–153):

| Кнопка | Клавиша | Кнопка | Клавиша |
|---|---|---|---|
| `16*f(a)` | `1` | `Next` | `N` |
| `8*f(a)` | `2` | `Back` | `B` |
| `8*f(psi)` | `3` | `Variables` | `V` |
| `2*f(a,t)` | `4` | `Constants` | `C` |
| `2*f(R,t)` | `5` | `Grids` | `D` |
| `8*f(t)` | `6` | `Save log` | `I` |
| `Equil` | `8` | `Land PS` | `Q` |
| `Layout` | `M` | `Port PS` | `G` |
| `Get X-axis` | `X` | `Write data` | `F` |
| `Style` | `.` | `Type model` | `L` |
| `Run` | Enter | `Step` | пробел |
| `Quit` | `/` | `Help` | `H` |

### 3.2. Обработчик `if_key`

Полный список ветвей `CASE` в `src/graph/gui_interaction.f90` (≈301–577):

| Клавиша | Действие |
|---|---|
| Enter (13) | продолжить расчёт |
| пробел (32) | пауза / один шаг (`TASK='DSP'`) |
| `%` | `cpu_report` — сводка затрат времени |
| `.` | стиль кривых |
| `/` | выход |
| `0`–`9` | режим графики; повторное нажатие `2`–`5` переключает `MODEY` (вариант раскладки), `6` перебирает 2/4 окна, `1` возвращает на первую вкладку |
| `H`, `?` | справка в терминал (28 строк `HELP`) |
| `N`, `B` | следующая/предыдущая вкладка |
| `C` | диалог констант (`CF`, `CV`, `CNB`, …) |
| `D` | диалог управляющих параметров и расписания подпрограмм (`DTEQ` и т. п.); `XOUT` 0–3 меняет ось X |
| `V` | диалог переменных; при изменении выставляет `IFDFVX=3` |
| `F` | `writeData` → `dat/<exp>-<equ>-<k>.dat` |
| `G` / `Q` | PostScript (портрет / ландшафт) → `dat/<exp>-<equ>-<k>.ps` |
| `I` | сохранить `equ/log/<equ>`: переменные до `ZRD1`, константы, 22 управляющих параметра |
| `J` | выполнить `ipcs -s` и `ipcs -m` (диагностика IPC, см. [07](07-heating-and-external-modules.md)) |
| `L` | напечатать `src/tmp/model.txt` (разобранная модель) |
| `M` | раскладка окон (`ASXWIN`, `ASTWIN`, `MENUTABLE` режима 7) |
| `W` | переразложить окна |
| `X` | информация о текущей оси X в терминал |
| Ctrl+буква | вызов подпрограммы модели, привязанной к этой клавише (`kibm2key`, `DTEQ(4,j)`) |
| Ctrl+C, Alt+/ | выход |
| Esc (27) | **не** выход: печатает `Use key "/" for exit` |

Ось X (`XOUT`, `writeData`): 0 — `a` на `[0, AB]`; 1 — `a` на `[0, ABC]`; 2 — `ρ` на `[0, ROC]`; 3 — `ψ` на `[FP(1), FP(NA1)]`.

### 3.3. Чего нет

Клавиши A7 `A` (цвета), `R`, `S`, `Y`, `T`, `U`, `P` в `if_key` A8 не обрабатываются, хотя текст справки всё ещё упоминает `A`, `Y`, `P`.

---

## 4. Режимы графики

`curves_per_frame = (16, 8, 8, 2, 2, 8, 4, 0, 0)` (`graph_utils.f90:158`); раскладка холста — `set_plot_area` (≈1605 и далее).

| Режим | Кнопка | Что показывает |
|---|---|---|
| 1 | `16*f(a)` | 16 радиальных профилей, холст 4×2 |
| 2 | `8*f(a)` | 8 профилей крупнее; повтор `2` — другая раскладка |
| 3 | `8*f(psi)` | 8 профилей по ψ |
| 4 | `2*f(a,t)` | 2 профиля с наложением по времени |
| 5 | `2*f(R,t)` | 2 профиля по большому радиусу R во времени |
| 6 | `8*f(t)` | временные ряды; 2 или 4 окна (`MODEY`) |
| 7 | — | используется раскладкой (`M`); 4 кривые |
| 8 | `Equil` | равновесие: LCFS, магнитные поверхности, стенка (`plot_lcfs`, `plot_flux_surfaces`, `plot_wall`) |
| 9, 0 | — | кривых нет (`curves_per_frame=0`) |

Пользовательской графической подпрограммы (A7 `mydraw`, режим 9) в A8 нет.

---

## 5. Экранный вывод в модели

### 5.1. Синтаксис

Разделы модели после уравнений (пример — `repo_user_area/equ/flux_feqis`):

```
	Profile output

Tex\TE\\TEX;	pinb\PIBM\-3;	Pecr\PEECR;	neut\NRATE;
q\1./MU\\CAR2X;	Prad\PRAD\-1;	...

	Time output

QE_QETOTB;	QI_QITOTB;	q0_1/MUC;	qa_1/MUB;
PNBi_VINT(PIBMB);	Pecr_VINT(PEECRB);	...
```

- Радиальный: `Имя\Выражение[\Масштаб]` или `Имя\Выражение\\Xdata[\Масштаб]`; `\\Xdata` добавляет экспериментальную кривую (точки) поверх расчётной; отрицательный масштаб объединяет окна в группы. Имя — до 4 символов. Не больше 128 строк (`config.NRW`; лишние отрезаются с предупреждением в `pyparse/equ_parser.py` ≈241–260).
- Временной: `Имя_Выражение[_Масштаб]` (`equ_parser.py` ≈218–231).
- Описание в документации A7 — `docu/section5.tex` ≈875–1072.

### 5.2. Куда это попадает

| Канал | Попадает ли экранный вывод |
|---|---|
| X11-окно | да |
| `F` → `dat/*.dat` | да, на текущий момент: временные выражения (`NAMET`/`TOUT`) и радиальные (`NAMER`/`ROUT`) |
| JSON / CDF | **нет** |
| batch-сборка | **не компилируется** вовсе |

Проверено на реальном CDF: выражения `q0`, `qa`, `QE_`, `PNBi`, `QETOTB` и т. п. в файле отсутствуют. Чтобы величина модели попала в CDF, её надо присвоить рабочему массиву или константе, которые ядро пишет всегда: `CARn` (профиль), `CFn` (скаляр), как это сделано в `flux_feqis`: `CAR1 = QETOT; CF15 = BETANB; CF8 = NEC/NEAVB;`. В CDF они появятся как `CAR1`, `CF15` и т. д. без единиц и описания.

---

## 6. JSON-вывод ядра

### 6.1. Когда пишется

`src/astra/astra_main.F90`: перед циклом `t_prev_json = TSTART - DPOUT`; на каждом проходе цикла по времени, если `TIME - t_prev_json > DPOUT`, вызывается `write_ajson`, затем `STEPUP`. `DPOUT` — управляющий параметр (по умолчанию 0.01 с, `src/astra/scalars.f90:128`; в `astra_variables.json` — «Profiles/time output interval»); в `equ/log/flux_feqis` задано `DPOUT=0.2`, поэтому для `-s 4 -e 5` получилось 6 снимков: 4.0, 4.2, …, 5.0.

### 6.2. Файл

`write_ajson` (`src/astra/json_rw.f90` ≈387–544) пишет `AWD/ncdf_out/<exp><equ>-<j_out>.json`, где `j_out = j_call + restart`. Размер одного файла для `flux_feqis` (NA1=91) — около 1.85 МБ.

Группы и число ключей (проверено на `aug34954flux_feqis-1.json`):

| Группа | Ключей | Форма значения | Содержимое |
|---|---|---|---|
| `variables` | 135 | скаляр | `AB`, `ABC`, `BTOR`, `IPL`, `ZRD1..`, … |
| `variables_x` | 135 | скаляр | экспериментальные аналоги `ABX`, `IPLX`, `ZRD1X`, … |
| `constants` | 162 | скаляр | `CF1..`, `CV1..`, `CNB1..`, … |
| `control` | 15 | скаляр | `DPOUT`, `TIME`, … |
| `internInt` | 33 | скаляр | `NEQUIL`, `MEQUIL`, `NA1`, … |
| `internDbl` | 46 | скаляр | `SGNIP`, `SGNBT`, … |
| `profiles_x` | 106 | `(NA1)` | `CAR1X`, `TEX`, `NEX`, … |
| `profiles` | 454 | `(NA1)` | `TE`, `NE`, `CU`, `MU`, `CAR1..`, `XRHO`, … |
| `equil_signals` | 8 | скаляр | `b0`, `betpol`, `i_plasma`, `li3`, `psibound`, `psiaxis`, `r0`, `Vloop` |
| `equil_profiles` | 28 | `(n_eq)` | `areat`, `bdb0`, `q`, `pprime`, `rho_tor_norm`, … |
| `equil_rz2d` | 5 | `(n_eq, n_theta)` | `r`, `rmin`, `psirz`, `theta2d`, `z` |
| `equil_coord` | 6 | `(n_eq, n_theta)` | `gradvcell`, `bpcell`, `bcell`, `rcell`, `darea`, `jphi` |
| `equil_rect` | 4 | `(nR, nZ)` | `psirz2d`, `fdia2d`, `r2d`, `z2d` |

Итого 1137 ключей — ровно столько переменных в CDF.

Особенности записи:

- имена и порядок берутся из `astra_variables.json` через `src/astra/json_vars.f90` (имена `character(len=6)`); `write_ajson` идёт по счётчику `jid`, поэтому **порядок записи в Fortran обязан совпадать с порядком в `astra_variables.json`** — при добавлении переменной надо править оба места;
- `write_array` обнуляет значения с `|x| < 1e-20` и `|x| > 1e20`;
- группа `strahl` (13 переменных) описана в `astra_variables.json`, но не пишется.

---

## 7. `json2cdf.py`: склейка в NetCDF

`postProcess/json2cdf.py`, функция `json_concat(expequ)`; вызывается из `repo_user_area/exe/as_exe.py` после `exe/Build` и очистки IPC, в `try: ... except: pass`. **Любая ошибка склейки молча проглатывается** — если CDF не появился или устарел, запускать вручную:

```bash
cd $AWD
python3 postProcess/json2cdf.py -e aug34954flux_feqis
```

Алгоритм:

1. читает `ncdf_out/<expequ>-1.json`, `-2.json`, … до первого отсутствующего файла **или** до файла, который старше (по `mtime`) предыдущего — так отбрасываются хвосты от прошлых, более длинных расчётов;
2. единицы и `long_name` берёт из `astra_variables.json` (поля `units`, `desc`);
3. определяет размерности: `TIME`, `XRHO` (=`NA1`), `RHO_SURF` (=`n_eq`), `THETA`, `R`, `Z`;
4. назначает размерности по форме массива и группе (`equil_*` → `RHO_SURF`/`THETA`/`R,Z`, остальные → `XRHO`);
5. пишет `ncdf_out/<expequ>.CDF` через `scipy.io.netcdf_file` (NetCDF3 classic), все данные `'>f8'`, включая целочисленные `internInt`.

Если в будущем появится переменная с формой, не попадающей ни в один из случаев шага 4, у неё не будет ключа `dimensions` и скрипт упадёт с `KeyError` (молча, см. выше).

**Дефект координатных переменных.** `RHO_SURF`, `THETA`, `R`, `Z` заполняются данными всех снимков сразу (`np.array(ds_astra['rho_tor_norm']['data'])` формы `(nt, n)`), хотя объявлены одномерными длины `n`. В файл записывается `nt·n` значений. Последствия, проверенные на реальном файле:

| Читатель | `RHO_SURF` | Почему |
|---|---|---|
| `netCDF4.Dataset` | `(91,)` — корректно | читает по длине размерности |
| `scipy.io.netcdf_file` | `(546,)` = 6×91 | берёт длину из записанного `vsize` |

Значения первых `n` элементов корректны (это сетка первого снимка; на проверенном файле сетка одинакова во всех снимках). При чтении через scipy надо обрезать: `v['RHO_SURF'].data[:f.dimensions['RHO_SURF']]`. `XRHO` и `TIME` записаны правильно.

---

## 8. Формат CDF

### 8.1. Имя и место

`ncdf_out/<exp><equ>.CDF` — конкатенация имени экспериментального файла и модели без разделителя: `aug34954` + `flux_feqis` → `aug34954flux_feqis.CDF`; `AUG33040_2500` + `fbe` → `AUG33040_2500fbe.CDF`. Каждый новый запуск той же пары перезаписывает файл.

### 8.2. Структура (реальный файл)

`~/astra/a8/ncdf_out/aug34954flux_feqis.CDF`: формат `NETCDF3_CLASSIC` (магия `CDF\x01`), ≈7 МБ.

| Размерность | Длина | Смысл |
|---|---|---|
| `TIME` | 6 | 4.0 … 5.0 с шагом `DPOUT=0.2` |
| `XRHO` | 91 | транспортная сетка, `NA1` |
| `RHO_SURF` | 91 | сетка магнитных поверхностей равновесия |
| `THETA` | 90 | полоидальный угол |
| `R`, `Z` | 65, 65 | прямоугольная сетка равновесия |

| Размерности переменной | Число | Примеры |
|---|---|---|
| `(TIME, XRHO)` | 559 | `TE`, `NE`, `TI`, `CU`, `MU`, `CAR1`, `TEX`, `PEBM` |
| `(TIME,)` | 534 | `IPL`, `BTOR`, `ZRD1`, `ZRD1X`, `CF1`, `DPOUT`, `NA1`, `psiaxis` |
| `(TIME, RHO_SURF)` | 27 | `q`, `pprime`, … |
| `(TIME, RHO_SURF, THETA)` | 10 | `r`, `z`, `rmin`, `psirz`, `gradvcell`, … |
| `(TIME, R, Z)` | 2 | `psirz2d`, `fdia2d` |

Плюс 5 одномерных координат `XRHO`, `RHO_SURF`, `THETA`, `R`, `Z` (сама `TIME` входит в 534); всего 1137. Атрибуты у каждой переменной — `units` и `long_name`.

| Переменная | `units` | `long_name` |
|---|---|---|
| `TE` | `keV` | Electron temperature |
| `NE` | `10^19/m^3` | Electron density |
| `CU` | `MA/m^2` | Longitudinal current density |
| `MU` | — | Inversed safety factor (1/q) |
| `q` | — | Safety factor from EQUIL |
| `QE` | `MW` | Total electron energy flux |
| `psiaxis` | `Wb` | Axis psi value |
| `r` | `m` | Major radius of magnetic surfaces' points |
| `psirz2d` | `Wb/rad` | Psi on rectangular grid, ASTRA has -2*PI*PSI_EQUIL |
| `BTOR`, `IPL` | `T`, `MA` | — |
| `TEX` | `keV` | External electron temperature |
| `CAR1`, `CF1` | пусто | рабочие массивы без описания |
| `XRHO` | пусто | Normalized Rho = Rho/Roc = [drhox/2,1-drhox/2] |

Значения на реальном файле: `XRHO[0]=0.00552`, `XRHO[-1]=1.00000004` — последняя точка равна 1, а не `1-drhox/2`, как утверждает `long_name`. `TE` при `t=4.6`: 4.087 кэВ в центре, 0.0764 кэВ на краю; `q0≈0.66`, `qa≈6.58` (`1/MU` на `XRHO`; `q` равновесия на `RHO_SURF` даёт 0.70 и 6.38 — сетки и определения границы разные); `IPL=0.797` МА; LCFS: R 1.15…2.13 м, Z −0.91…0.75 м.

### 8.3. Связь с MDSplus/IMAS

Готового экспорта в MDSplus в дереве нет. Есть `postProcess/acdf2imas.py` (CDF → IMAS, §11). Для собственного импорта CDF удобен тем, что все переменные однородны (`>f8`), с единицами, на 5 фиксированных размерностях.

---

## 9. Чтение CDF на Python

`requirements.txt` гарантирует только `scipy`, поэтому основной пример — на `scipy.io.netcdf_file`. Пример проверен на реальном файле (scipy 1.11.4, numpy 1.26.4):

```python
import numpy as np
from scipy.io import netcdf_file

fcdf = 'ncdf_out/aug34954flux_feqis.CDF'
with netcdf_file(fcdf, 'r', mmap=False) as f:
    v = f.variables
    t   = v['TIME'].data.copy()             # (nt,) с
    rho = v['XRHO'].data.copy()             # (91,) ρ/ROC
    te  = v['TE'].data.copy()               # (nt, 91) кэВ
    print(v['TE'].units.decode(), v['TE'].long_name.decode())   # keV Electron temperature

    jt = np.argmin(np.abs(t - 4.6))
    q = 1.0 / v['MU'].data[jt]              # q на транспортной сетке
    ipl = v['IPL'].data.copy()              # (nt,) МА

    n_eq = f.dimensions['RHO_SURF']         # обрезка из-за дефекта json2cdf (§7)
    rho_surf = v['RHO_SURF'].data[:n_eq].copy()
    q_eq = v['q'].data[jt]                  # q из равновесия, по RHO_SURF

    rb = v['r'].data[jt, -1, :]             # граница плазмы (последняя поверхность)
    zb = v['z'].data[jt, -1, :]
    R = v['R'].data[:f.dimensions['R']].copy()
    Z = v['Z'].data[:f.dimensions['Z']].copy()
    psi = v['psirz2d'].data[jt]             # (65, 65), оси (R, Z)
```

Замечания:

- атрибуты в scipy — `bytes`, отсюда `.decode()`;
- с `mmap=False` массивы уже в памяти, но `.copy()` полезен, если файл закрывается до использования данных;
- профили (`XRHO`) и величины равновесия (`RHO_SURF`) лежат на **разных** сетках; совпадение длин (91 и 91) — случайность конкретного запуска, их нельзя смешивать без интерполяции;
- `psirz2d` индексируется `[время, R, Z]`.

Альтернатива — `netCDF4` (если установлен): `netCDF4.Dataset(fcdf)['TE'][:]` возвращает `numpy.ma.MaskedArray`; координаты читаются правильно без обрезки. `xarray` в окружении машины разработчика не установлен; с ним `xarray.open_dataset(fcdf)` должен работать через бэкенд `netcdf4` или `scipy` (не проверено).

Перечень переменных по размерностям:

```python
from collections import defaultdict
groups = defaultdict(list)
with netcdf_file(fcdf, 'r', mmap=False) as f:
    for name, var in f.variables.items():
        groups[var.dimensions].append(name)
for dims, names in groups.items():
    print(dims, len(names))
```

---

## 10. Рестарт из JSON

`as_exe ... -re N` записывает `restart=N` в `tmp/astra.nml`; `astra_main.F90` вместо начальных итераций вызывает `read_ajson(N)` — состояние берётся из `ncdf_out/<exp><equ>-<N>.json`, `TEND` — из командной строки. Нумерация новых JSON продолжается с `N+1` (`j_out = j_call + restart`), а `json2cdf.py` по-прежнему склеивает с `-1.json`: старые файлы `1..N` старше новых, поэтому проверка `mtime` их пропускает. Это заменяет механизм A7 `RESTARTA/IRESTA`. Рестарт в GUI-режиме — также единственный способ «посмотреть» на сохранённый момент расчёта в окне ASTRA.

---

## 11. Скрипты `postProcess/` и `plot_compare.py`

| Скрипт | Назначение | Зависимости | Замечания |
|---|---|---|---|
| `json2cdf.py` | JSON → CDF | scipy | §7 |
| `acdf2imas.py` | CDF → IMAS IDS `core_profiles`, `equilibrium`, `core_sources` | модуль `imas`, `IMASDB`; бэкенды HDF5/MDS+/ASCII | файл по умолчанию `../ncdf_out/aug34954fluxes.CDF` — модели `fluxes` в дереве нет |
| `eqdsk.py` | класс `EQDSK`: чтение/запись/графики, COCOS (`to_coco`, `find_coco`), `Brzt` | `rw_for` (рядом), matplotlib | лимитер AUG зашит; `__main__` ссылается на путь IPP |
| `json2eqdsk.py` | JSON-снимок → G-EQDSK (nR=129, nZ=257, `cocos_out=7`; Delaunay + CloughTocher) | scipy, `eqdsk.py` | **дефект**: обращается к `json_d['internal']['SGNIP']`, а группа называется `internDbl` → `KeyError: 'internal'` (проверено запуском на реальном JSON). Альтернатива — `EQDSK(coco)` в модели (§1) |
| `plotdkesrun.py` | разбор вывода DKES, `output.txt` + PNG | pandas, matplotlib | стелларатор |
| `readxmxn.py` | чтение мод `xm`/`xn` из `dat/wout_VMECoutput.nc` | netCDF | стелларатор |
| `rw_for.py` | форматированный ввод/вывод Fortran (`wr_for`, `ssplit`, `lines2fltarr`, `fltarr_len`) | numpy | копия есть в `pyparse/rw_for.py` |
| `test_pressure.py` | сравнение `NE*TE+NI*TI+`давление быстрых с `pressure` равновесия | scipy, matplotlib | путь `../ncdf_out/aug34954flux_feqis.CDF` относительный — запускать из `postProcess/` |

`plot_compare.py` (корень): рисует `pprime` в последний момент из `regressions/<f>.CDF` и `ncdf_out/<f>.CDF`, где `fname='aug34954fluxes'` зашит — такого эталона нет, без правки скрипт падает. Кроме того, `pprime` лежит на `RHO_SURF`, а рисуется против `XRHO`; это работает только потому, что обе длины 91.

При импорте этих скриптов из дерева исходников стоит использовать `python3 -B` (или `PYTHONDONTWRITEBYTECODE=1`), чтобы не создавать `__pycache__` внутри `astra-src`.

---

## 12. Регрессии: `regressions/`, `compReg.sh`, `compareRegressions.py`

### 12.1. Эталоны

Имя эталона — то же `<exp><equ>.CDF`, что и у расчёта.

| Файл в `regressions/` | `exp` | `equ` | Особенности |
|---|---|---|---|
| `AUG33040_2500fbe.CDF` | `AUG33040_2500` | `fbe` | 15 моментов 2.48–2.62; без нагрева |
| `AUG36982_3400tglf_pid.CDF` | `AUG36982_3400` | `tglf_pid` | `XRHO`=161, t 4.0–4.7; `PEECR` до 10.9 (TORBEAM) |
| `aug34954_tflux_feqis.CDF` | `aug34954_t` | `flux_feqis` | `TE` из неусреднённого `E34954.IDA`, без `ZRD3`/`ZRD8` |
| `aug34954astrahl_simple_msp.CDF` | `aug34954` | `astrahl_simple_msp` | **устаревший** (коммит `7a4d9a19`, 2026-07-14), 1122 переменных; нет `FTPT`, `FTPTX`, `GVAC`, `GVACX`, `PHIEDG`, `PHIEDGX`, `SG11..SG22(X)`, `vmec_option`; в `compReg.sh` не входит |
| `aug34954flux_cuas.CDF` | `aug34954` | `flux_cuas` | — |
| `aug34954flux_feqis.CDF` | `aug34954` | `flux_feqis` | `PEECR` 5.94, `PIBM` 1.34, `PIICR` 0.28 — посчитан с TORBEAM и RABBIT |
| `aug34954flux_nbi.CDF` | `aug34954` | `flux_nbi` | `PIBM` 0.76 (встроенный NBI), `PEECR` от TORBEAM |
| `aug34954flux_neo_tglf.CDF` | `aug34954` | `flux_neo_tglf` | 21 момент; без нагрева |
| `aug34954flux_spider.CDF` | `aug34954` | `flux_spider` | SPIDER |
| `aug34954qlk.CDF`, `aug34954qlknn.CDF`, `aug34954tglf.CDF` | `aug34954` | `qlk`, `qlknn`, `tglf` | `tglf` — 21 момент, 23.9 МБ |

Все, кроме `astrahl_simple_msp`, обновлены в коммите `e66f6ba1` (2026-10-01).

### 12.2. `compReg.sh`

Подключает `platform/env.<pl>` и для каждой пары вызывает `$AWD/exe/as_exe -m <equ> -v <exp> -s .. -e ..`, затем `python3 $AWD/compareRegressions.py -m <equ> -v <exp>`:

| `exp` | Модели | Интервал |
|---|---|---|
| `aug34954` | `flux_spider`, `flux_cuas`, `flux_nbi`, `flux_feqis`, `qlk`, `tglf`, `qlknn` | 4–5 с |
| `aug34954_t` | `flux_feqis` | 4–5 с |
| `AUG33040_2500` | `fbe` | 2.48–2.7 с |
| `AUG36982_3400` | `tglf_pid` | 4–6 с |
| `aug34954` («slow») | `flux_neo_tglf` | 4–5 с |

Замечания по коду:

- `as_exe` вызывается **без `-batch`**, т. е. в GUI-режиме: нужен X-дисплей и GUI-сборка (окно само закроется в конце);
- `compareRegressions.py` открывает `ncdf_out/...` и `regressions/...` по **относительным** путям, поэтому `compReg.sh` надо запускать из `$AWD`;
- скрипт не возвращает код ошибки при расхождениях — результат только в логах;
- на машине без закрытых модулей и с `MAX_NWORKERS<64` большая часть моделей остановится (заглушки, IPC), см. [07](07-heating-and-external-modules.md).

### 12.3. `compareRegressions.py`

`compare(fcdf, fcdf_ref, tolerance=1e-7)`:

1. сообщает о ключах эталона, которых нет в новом файле, и наоборот;
2. для общих ключей при несовпадении формы — `Shape mismatch`;
3. иначе для каждого момента `jt`:

\[
\mathrm{diff} = \frac{\lVert x_\text{new} - x_\text{ref} \rVert_2}{\sum |x_\text{new}|}
\]

(если `sum|x_new|=0` — нормировка на `sum|x_ref|`; если обе нулевые — 0). При `diff > tolerance` пишется `ERROR: <ключ>: discrepancy ... at time=...`.

Числитель — евклидова норма, знаменатель — L1-норма, поэтому для профиля из `n` точек `diff` примерно в `√n` раз меньше «настоящей» относительной ошибки: допуск `1e-7` для профиля на 91 точке соответствует относительной ошибке порядка `1e-6`.

Лог: `regressions/<exp><equ>.log` (каталог скрипта), уровень `ERROR`, перезаписывается. Без аргументов по умолчанию `-m fluxes -v aug34954` — эталона `aug34954fluxes.CDF` нет.

### 12.4. Известные расхождения

Из `astra-src/AGENTS.md` и опыта сборки на Ubuntu (`installer/docs/DETAILED.md`):

- gfortran против Intel + MKL: допуск `1e-7` **всегда** превышается — эталоны посчитаны компилятором Intel с MKL;
- `fbe`/`AUG33040_2500`: расхождения ≈0.1 % по `CU`/`MU`, 1–3 % по `psiaxis`; на последнем шаге `t=2.62` — `Error in find new boundary`, так же, как в эталоне;
- `flux_feqis`/`aug34954` без RABBIT и TORBEAM не сравним с эталоном (нет `PEECR`, `PIBM`, `PIICR` от закрытых кодов);
- OpenBLAS против MKL — расхождения до `3e-6`.

Практически на gfortran сравнение осмысленно делать с собственным эталоном (`regressions/` своей сборки) или с увеличенным допуском, вызывая `compare(..., tolerance=1e-3)` из Python.

---

## 13. Документация A7 против реальности A8

`docu/section5.tex` ≈2814–3440 («Brief operation guide») описывает ASTRA 7.

| Тема | A7 (`docu/section5.tex`) | A8 |
|---|---|---|
| Запуск | `Astra [data] [model] [t0] [t1] [postview] [test]`, опции `-v -m -s -e -r` | `exe/as_exe -v -m -s -e [-batch] [-re N] ...` ([02](02-build-and-run.md)) |
| Просмотр после расчёта | `Review`, режим postview | нет; ближайшая замена — рестарт `-re N` в GUI или чтение CDF |
| Клавиши | в т. ч. `A`, `R`, `S`, `Y`, `T`, `U`, `P`, Esc | нет; выход только `/` |
| Пользовательская графика | `mydraw`, режим 9 | нет |
| Вывод | U-файлы, `dat/` | JSON → CDF; U-файлы не пишутся |
| Синтаксис экранного вывода | `Name\Expression\\Xdata\Scale`, `Name_Expression_Scale` | тот же, обрабатывает `pyparse/equ_parser.py` |

---

## 14. Известные дефекты и подводные камни

| Что | Где | Последствие |
|---|---|---|
| Ошибки `json2cdf` проглатываются | `exe/as_exe.py`, `try/except: pass` | CDF может не обновиться без сообщений |
| Координаты `RHO_SURF`/`THETA`/`R`/`Z` записаны с данными всех моментов | `postProcess/json2cdf.py` | scipy читает `nt·n` значений; обрезать по размерности |
| `long_name` у `XRHO` не соответствует данным | `astra_variables.json` | последняя точка — 1.0 |
| Экранный вывод не сохраняется | архитектура (`astra_out.f90` только в GUI) | присваивать в `CARn`/`CFn` |
| `json2eqdsk.py`: `json_d['internal']` | `postProcess/json2eqdsk.py` | `KeyError`; использовать `EQDSK(coco)` |
| Устаревший эталон `astrahl_simple_msp` | `regressions/` | при сравнении будут `Key ... is missing` |
| `plot_compare.py`, `compareRegressions.py`, `acdf2imas.py`: имя `aug34954fluxes` | корень, `postProcess/` | нет такой модели/эталона |
| Порядок переменных JSON | `json_rw.f90` против `astra_variables.json` | при рассинхронизации — сдвиг имён |
| `as_exe` удаляет только сегменты памяти IPC, семафоры остаются | `exe/as_exe.py` | `ipcs -s` / `ipcrm -s` вручную |
| `exe/Build` печатает `ASTRA normal exit` безусловно | `exe/Build` | признак успеха — строка `>>> ASTRA normal exit >>>` от ядра |

---

## 15. Где смотреть в исходниках

| Что | Где |
|---|---|
| Главный цикл, вызовы `gui_init`, `write_ajson`, `read_ajson` | `src/astra/astra_main.F90` |
| Запись/чтение JSON | `src/astra/json_rw.f90` (`write_ajson` ≈387–544, `read_ajson`) |
| Имена и группы переменных | `astra_variables.json`, `src/astra/json_vars.f90` |
| `DPOUT` по умолчанию | `src/astra/scalars.f90:128` |
| Склейка в CDF | `postProcess/json2cdf.py` |
| Вызов склейки и очистка IPC | `repo_user_area/exe/as_exe.py` |
| X11-сборка | `repo_user_area/exe/Makefile` ≈52–72 |
| Генерация `RADOUT`/`TIMOUT`/имён | `pyparse/code_gen.py` ≈311–358, `pyparse/const_text.py` ≈586–624, `pyparse/equ_parser.py` ≈218–260 |
| Низкоуровневая графика | `src/graph/Astra2XW.c` |
| Меню и диалоги | `src/graph/dialogx.c` (кнопки ≈140–153) |
| Режимы, окна, раскладка | `src/graph/graph_utils.f90` (`gui_init`, `plotMode`, `set_plot_area`, `curves_per_frame` строка 158) |
| Клавиши, `writeData`, PS, лог | `src/graph/gui_interaction.f90` (`if_key` ≈301–577, `writeData` ≈714, `set_filename` ≈1334) |
| EQDSK из модели | `repo_user_area/sbr/a2eqdsk.f90` |
| Скрипты постобработки | `postProcess/*.py`, `plot_compare.py` |
| Регрессии | `regressions/*.CDF`, `compReg.sh`, `compareRegressions.py` |
| Известные расхождения | `AGENTS.md` (корень `astra-src`), `installer/docs/DETAILED.md` |
| Документация A7 | `docu/section5.tex` ≈875–1072 (экранный вывод), ≈2814–3440 (руководство оператора) |
