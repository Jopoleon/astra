# 09. Каталог пользовательской области (`repo_user_area/`)

Состояние: исходники ASTRA 8 в `astra-src/`, ветка `local`, коммит `67281112` (2026-10-01, «ASTRA 8.6, January 2026» по `exe/version`).

Глава описывает всё, что лежит в `repo_user_area/`: модели (`equ/`), эксперименты (`exp/`), данные (`udb/`), формулы (`fml/`), функции (`fnc/`), подпрограммы (`sbr/`), а также `xpr/`, `preProcess/`, `pyparse/` и `exe/`. Описания составлены по содержимому файлов (заголовки, комментарии, код) и по `docu/section4.tex` (разделы «Astra expressions → List of expressions», «Built-in functions», «Plug-in subroutines»). Где назначение по файлу определить нельзя, так и написано: «назначение не установлено».

Связанные главы:

- синтаксис моделей (`:EQ`, `:AS`, `<:`, `>::::`, вывод `\` и `_`) — [03-model-language.md](03-model-language.md);
- формат файлов эксперимента и U-файлов — [04-experiment-data.md](04-experiment-data.md);
- решатель уравнений переноса — [05-solver-core.md](05-solver-core.md);
- равновесие (`IPEQL`, FEQIS, SPIDER, EMEQ, VMEC) — [06-equilibrium.md](06-equilibrium.md);
- нагрев и внешние модули (RABBIT, TORBEAM, NBI, TORIC, TGLF, NEO, QuaLiKiZ, STRAHL) — [07-heating-and-external-modules.md](07-heating-and-external-modules.md);
- графика и вывод (`Profile output`, `Time output`, CDF) — [08-graphics-and-output.md](08-graphics-and-output.md);
- сборка и запуск (`exe/`, `as_exe`, `Build`, `astra_rc`) — [02-build-and-run.md](02-build-and-run.md);
- термины — [glossary.md](glossary.md).

---

## 1. Как пользовательская область попадает в рабочий каталог

`install.sh` в корне исходников копирует каждый подкаталог `repo_user_area/<DIR>` в `$AWD/<DIR>`:

- без ключа `-safe` — копирует всё молча (перезаписывая);
- с ключом `-safe` — `equ/`, `exp/`, `udb/` и отсутствующие каталоги копируются всегда, для остальных спрашивается «backup + copy»; дальше задаются вопросы про RABBIT, TORBEAM, QuaLiKiZ, QLKNN, NEO, TGLF. Ответ `n` про RABBIT/TORBEAM комментирует `RABBIT`/`TORBA` в `equ/flux_feqis` (`sed`), ответ `n` про QLK/NEO/TGLF убирает цели из строки `all:` в `exe/Makexpr`.

В конце `install.sh` всегда делает `make -f exe/Makefile clean` (без `-safe`) и запускает тестовый расчёт `exe/as_exe -m flux_feqis -v aug34954 -s 4 -e 5`.

Наш установщик (`installer/install_astra8.sh`, монорепозиторий) делает то же копирование, но дополнительно правит рабочую копию: комментирует `^TORBA`/`^RABBIT` только в `equ/flux_feqis` и убирает TGLF/QLK/NEO из `exe/Makexpr` (см. [02-build-and-run.md](02-build-and-run.md)).

Состав `repo_user_area/`:

| каталог | файлов | содержание |
|---|---|---|
| `equ/` | 27 моделей + 26 файлов `equ/log/` = 53 | модели и стартовые значения констант |
| `exp/` | 6 файлов эксперимента + `readme` + подкаталоги `cnf/` (4), `nbi/` (1), `nml/` (8), `icrh/` (1) | описания разрядов |
| `udb/` | 3 подкаталога, 32 U-файла | данные ASDEX Upgrade |
| `fml/` | 311 | формулы (фрагменты Fortran, вставляемые в генерируемый код) |
| `fnc/` | 117 | функции (`*.f90`, собираются в `lib/libfnc.a`) |
| `sbr/` | 37 | подпрограммы (собираются в `lib/libsbr.a`) |
| `xpr/` | 4 | IPC-обёртки TGLF / QuaLiKiZ / NEO |
| `preProcess/` | 2 | подготовка входных данных из IMAS / shotfiles AUG |
| `pyparse/` | 1 | дополнение к парсеру моделей |
| `exe/` | 14 | скрипты сборки и запуска |

---

## 2. Модели `equ/`

### 2.1. Общие сведения

Модель — текстовый файл на языке ASTRA (см. [03-model-language.md](03-model-language.md)). Для каждой модели `<m>` обязателен файл `equ/log/<m>`: `src/astra/read_input.f90` (`readInput`) читает его при старте и останавливается с `>>> Error: file "equ/log/<m>" missing`, если его нет. Для модели в подкаталоге `equ/<sub>/<m>` ищется `equ/<sub>/log/<m>`.

Выбор решателя равновесия задаётся переменной `IPEQL` (по умолчанию `5`, см. `src/astra/scalars.f90`; ветвление — `src/astra/metrics.f90`, подпрограмма `METRIC`):

| `IPEQL` | решатель | примечание |
|---|---|---|
| `-2` | цилиндрическая геометрия (`EQCYL`) | без тороидальности |
| `-1` | метрика из файла эксперимента (`set_external_metric`) | |
| `-6`, `6` | внешняя метрика (`set_external_metric_2`), `6` — для стелларатора | |
| `0` | `EQGUESS` (нет решателя, `NEQUIL=0`) | |
| `1` | EMEQ (`A2EMEQ`, `src/astra/emeq.f90`) | открытый, встроенный |
| `3` | только итерации `RHSEQ` | |
| `4` | SPIDER (`A2GSSOLVER`, `equil_solver=3`) | закрытый код; вызов `spider_run` в `src/astra/gs_solver.F90` обёрнут в `#ifdef SPIDER` |
| `5` | FEQIS (`A2GSSOLVER`, `equil_solver=101`) | открытый, `src/feqis/` |
| `7` | VMEC через STELLOPT (`vmec_interface`), режим задаёт `vmec_option` 0–9 | стелларатор; нужен `$STELLOPT_HOME` |

Все `sbr/*` компилируются в `lib/libsbr.a` целиком (`exe/Makesub`, wildcard по `*.f90/*.F90`). Для закрытых и внешних кодов в `sbr/` есть заглушки под `#ifdef`: `torba.F90` (`TORBEAM_INSTALLED`), `rabbit.F90` (`RABBIT_INSTALLED`), `torfpql_mod.F90` (`TORIC_INSTALLED`), `tglf_serial.F90` (`TGLF_INSTALLED`), `qlknn_serial.F90` (`QLKNN_INSTALLED`). Заглушка печатает `... not installed, check ..._LIB path in exe/astra_rc` и делает `stop`. Значит, модель с незакомментированным `TORBA`/`RABBIT` **собирается, но останавливается при первом вызове** — поэтому такие строки нужно комментировать `!`.

IPC-вызовы `tglf_ipc`, `neo_ipc`, `qlk_ipc` (`sbr/*_ipc.f90`) запускают отдельные программы `xpr/tglfi`, `xpr/neo`, `xpr/qlki`, которые собирает `exe/Makexpr` из `xpr/*` и библиотек TGLF/NEO/QuaLiKiZ (нужен MPI). Установщик эти цели убирает, поэтому IPC-модели на чистой Ubuntu не работают.

Сопоставление с эталонами `regressions/*.CDF`: имя эталона = `<эксперимент><модель>.CDF`.

| эталон | эксперимент | модель |
|---|---|---|
| `AUG33040_2500fbe.CDF` | `AUG33040_2500` | `fbe` |
| `AUG36982_3400tglf_pid.CDF` | `AUG36982_3400` | `tglf_pid` |
| `aug34954_tflux_feqis.CDF` | `aug34954_t` | `flux_feqis` |
| `aug34954astrahl_simple_msp.CDF` | `aug34954` | `astrahl_simple_msp` |
| `aug34954flux_cuas.CDF` | `aug34954` | `flux_cuas` |
| `aug34954flux_feqis.CDF` | `aug34954` | `flux_feqis` |
| `aug34954flux_nbi.CDF` | `aug34954` | `flux_nbi` |
| `aug34954flux_neo_tglf.CDF` | `aug34954` | `flux_neo_tglf` |
| `aug34954flux_spider.CDF` | `aug34954` | `flux_spider` |
| `aug34954qlk.CDF` | `aug34954` | `qlk` |
| `aug34954qlknn.CDF` | `aug34954` | `qlknn` |
| `aug34954tglf.CDF` | `aug34954` | `tglf` |

Эталоны получены на кластере IPP (Intel + MKL) с установленными RABBIT/TORBEAM/SPIDER/TGLF и т. д., поэтому для моделей с этими кодами расчёт на Ubuntu с эталоном не сравним.

### 2.2. Семейство `flux_*` (интерпретационные расчёты AUG #34954)

Модели `flux_*` построены по одному шаблону (базовый — `flux_feqis`). Общее для всех:

- профили `NE`, `TE`, `TI` **заданы** из эксперимента (`NE:AS; NE=NEX;` и т. д.), решается только **уравнение диффузии тока** (`CU:EQ`) или ток тоже задан (`CU:AS`);
- начальный `MU` (1/q) — из `CAR2X` (в `exp/aug34954` это U-файл `Q34954.QPR`, профиль q);
- ионный состав: дейтерий (`AMJ`, `ZMJ`) + углерод (`AIM1=12`, `ZIM1=ZICAR`, `NIZ1=CIMP1*NE`), `ZEF` из квазинейтральности;
- излучение: `PBOL1 = PRBER*NIZ1*NE` (комментарий в модели говорит «radiation due to C», но формула `fml/prber` — для бериллия, Z=4, A=9), плюс тормозное `PBRAD` и синхротронное `CRAD2*PSYNC`;
- нейтралы: `NEUT:.05:::;` (подпрограмма `sbr/neut.f90`, открытая);
- неоклассика: проводимость Sauter `CC=CNSA`, бутстреп по Kim (`HC=HCKIM; DC=DCKIM; XC=XCKIM`), источник тока `CD=CUBM+CUECR` (NBI + ECCD);
- полоидальная скорость `VPOL=VPSWW`, тороидальная `VTOR=VTORX`, `ER` из радиального баланса сил;
- сетка равновесия `NEQUIL=91`, `MEQUIL=90`, шаг по времени `TAUMIN=TAUMAX=0.05` с.

Ниже — только отличия.

**`flux_feqis`** — эталонная тестовая модель (запускается `install.sh` и нашим установщиком). `IPEQL=5` (FEQIS). Решается только `CU:EQ`. Вызывает `TORBA<:.05:::;` (TORBEAM → `PEECR`, `CUECR`) и `RABBIT(pRF_MW=1, nRF_harm=2, fRF_MHz=36.5, pi_icr=PIICR, pe_icr=PEICR)<:.05:::;` (RABBIT → `PEBM`, `PIBM`, `SNEBM`, `CUBM`, плюс простая модель ИЦР-нагрева с выводом в `PIICR`/`PEICR`). Добавлен водородный миноритет `NHYDR=0.05*NDEUT`. Эксперимент — `aug34954` (эталоны `aug34954flux_feqis.CDF` и `aug34954_tflux_feqis.CDF`). Нужны RABBIT и TORBEAM (закрытые IPP); после комментирования двух строк работает на Ubuntu (проверено установщиком, 30–45 с на `-s 4 -e 5`).

**`flux_cuas`** — «flux, CU assigned»: ток **задан** (`CU:AS`), строки `CC`/`CD` закомментированы, уравнение тока не решается. Дополнительно `EQDSK(7)<:0.05:::;` — запись файла G-EQDSK в соглашении COCO 7 (`sbr/a2eqdsk.f90`, открытый). Во временном выводе сравнение `Ipl1_IINT(CUB)` и `Ipl2_IPL`. `IPEQL=5`. Вызывает TORBA и RABBIT (без ИЦР-аргументов). Эксперимент `aug34954` (эталон есть).

**`flux_emeq`** — то же, что `flux_feqis`, но с `IPEQL=1` (встроенный трёхмоментный решатель EMEQ) и без расчёта `ER`. Вызывает TORBA и RABBIT. Эксперимент — по аналогии `aug34954` (эталона нет).

**`flux_nbi`** — вместо RABBIT используется открытый пакет NBI А. Р. Полевого (`sbr/nbi.f90` → `src/nbi/`): `CNB1=8; CNB2=1; CNB3=2; CNB4=1; NBI<:0.05:::;`, геометрия 8 инжекторов берётся из `exp/nbi/aug34954`. Дополнительно `TORBA<:0.05:::;` и `EQDSK(7)`. `NEUT<:0.05:::;` (с `<`). `IPEQL=5`. Эксперимент `aug34954` (эталон есть). Закрытый код только TORBEAM.

**`flux_neo`** — диагностический вызов NEO (GACODE) через IPC: `NEO_IPC<:0.1:::;`, результаты `neo_out%chi_i`, `chi_e`, `e_pflux` выводятся как `CAR17`, `CAR18`, `CAR20` (только на графики, в перенос не подставляются). TORBEAM/RABBIT закомментированы. `IPEQL=5`. Нужен NEO + MPI (`xpr/neo`).

**`flux_neo_tglf`** — то же, но вызываются и `TGLF_IPC<:0.1:::;`, и `NEO_IPC<:0.1:::;`; результаты в `CAR17–CAR18, CAR20` (TGLF) и `CAR21–CAR23` (NEO). `IPEQL=4` (SPIDER). Эксперимент `aug34954` (эталон есть). Нужны SPIDER, TGLF, NEO.

**`flux_pel`** — моделирование пеллет-инжекции: `ABLATION(CF16, CAR51)>::::;` (`sbr/ablation.f90`) даёт профиль абляции, источник `CAR53 = CF16*CAR51/VINT(CAR51B)` добавляется в `SN=SNEBM+CAR53`. Шаг по времени 5 мс. Вызывает TORBA и RABBIT. `ABLATION` читает namelist `&pellet` (`rho_abl_file`, `time_abl_file`, `mass`) из `exp/nml/<эксперимент>` и два U-файла — в репозитории ни в одном `exp/nml/*` группы `&pellet` нет, так что данных для запуска нет. Эксперимент не установлен.

**`flux_spider`** — `flux_feqis` с `IPEQL=4` (SPIDER) и `EQDSK(7)`. Эксперимент `aug34954` (эталон есть). Нужны SPIDER, TORBEAM, RABBIT.

**`flux_tglf`** — диагностический вызов `TGLF_IPC<:0.1:::;` на заданных профилях; `chi_i`, `chi_e`, `e_pflux` → `CAR17`, `CAR18`, `CAR20` (графики `chii`, `chie`). TORBEAM/RABBIT закомментированы. `IPEQL=5`. Нужен TGLF + MPI.

**`flux_tglf_serial`** — то же, но через `tglf_serial(CAR17, CAR18, CAR20, CAR42, CAR31, …, CAR35)>:0.1:::;` — TGLF, слинкованный прямо в исполняемый файл (`sbr/tglf_serial.F90`, `#ifdef TGLF_INSTALLED`), без IPC. Нужна библиотека TGLF.

**`flux_toric`** — `flux_feqis` плюс ИЦР-нагрев полноволновым кодом TORIC-SSFPQL: `TORIC(pwic=1.0d0, nmax=2, debug=0)<:.05:::;` (`sbr/torfpql_mod.F90`, обёртка R. Bilato). TORIC читает шаблон `exp/icrh/torica.inp` и группу `&toric` из `exp/nml/aug34954`. RABBIT вызывается с `pRF_MW=2` и выводом ИЦР в `CAR11`/`CAR12`. Нужны TORIC, RABBIT, TORBEAM. Путь к TORIC в `exe/astra_rc` указывает на дерево `/mpcdf/...` кластера IPP.

### 2.3. Предиктивные модели с TGLF / QuaLiKiZ

Общий шаблон (`tglf`, `qlk`, `qlknn`, `impli`, `tglf_*`): «AUG», `IPEQL=5`, `ATREQ=1`; **решаются** уравнения `NE`, `TE`, `TI` (`:EQ[2,CF16]`, т. е. до границы `ρ=CF16=0.85`, снаружи — эксперимент); транспорт: турбулентный χ из TGLF/QuaLiKiZ (ограничен `min(20, max(0, …))`) плюс неоклассика (`HNASI` для ионов, `HNGSE` для электронов), `DN=HNGSE`, `CN` = турбулентный поток частиц минус пинч Уэра `VP*VRAS`; бутстреп Sauter (`HCSA/DCSA/XCSA`), проводимость Hirshman `CC=CNHR`; `PRAD=PBRAD`; нагрев TORBEAM + RABBIT (кроме отмеченных). Граничные «добавки» к диффузии `DVN=DVE=DVI=CF14`.

**`tglf`** — `tglf_ipc(CF16)<:0.02:::;`, используются `chi_i`, `chi_e`, `e_pflux`, `mom_flux`, `equipart`, `gamma`, `omega`. `CU:EQ`. Эксперимент `aug34954` (эталон есть). Нужны TGLF (+MPI), TORBEAM, RABBIT.

**`tglf_pid`** — как `tglf`, но с обратной связью по плотности: `PID_CONTROL(…, CF1-CF2, NNCL, CF3)>::::;` (`sbr/pid_control.f90`) подстраивает плотность холодных нейтралов `NNCL`, чтобы средняя по объёму плотность `CF2=VINT(NEB)/VOLUME` следовала за экспериментальной `CF1`; затем `NEUT>::::;`. TGLF на всём радиусе (`CF16=1.0`). Эксперимент `AUG36982_3400` (эталон есть). Нужны TGLF, TORBEAM, RABBIT.

**`tglf_rot_pred_08`** — как `tglf`, плюс **предсказание вращения**: решается `UPAR:EQ[2, CF16]` с потоком импульса от TGLF (`CAR45`, сглаживание `SMEARR2`), числом Прандтля `XUPAR=20`, пинчем `CNPAR`; момент от NBI `TTRQ = RTOR*SCUBM*5.977e7` (комментарий: «rtor depends on RABBIT output»); `VTOR=UPAR`. Значение суффикса `_08` не установлено. Эксперимент не установлен (по комментарию «AUG»). Нужны TGLF, TORBEAM, RABBIT.

**`tglf_serial`** — как `tglf`, но через `tglf_serial(…)<:0.02:::;` без IPC; обмен энергией из TGLF (`CAR33`) добавлен в `PE` и вычтен из `PI`. Нужны библиотека TGLF, TORBEAM, RABBIT.

**`impli`** — как `tglf` (вызов `tglf_ipc<:0.02:::;` без аргумента), но уравнения `TE`/`TI` записаны со звёздочкой (`TE*:EQ`, `TI*:EQ`) — по `pyparse/eqns_init.py` это включает совместное неявное решение `TE`–`TI` (`eqns.tetieqn`); граничные значения `NEB`, `TEB`, `TIB` берутся из эксперимента в точке `ρ=CF16`. `IPEQL` не задан (по умолчанию 5, FEQIS). Эксперимент не установлен. Нужны TGLF, TORBEAM, RABBIT.

**`qlk`** — как `tglf`, но с QuaLiKiZ: `qlk_ipc(CF16)<:0.02:::;`, χ сглаживаются `SMEARR2`; ток задан (`CU:AS`). Эксперимент `aug34954` (эталон есть). Нужны QuaLiKiZ (+MPI), TORBEAM, RABBIT.

**`qlknn`** — как `qlk`, но с нейросетевой аппроксимацией QLKNN: `qlknn_serial(CAR21, CAR22, CAR24)<::::;` (`sbr/qlknn_serial.F90`, `#ifdef QLKNN_INSTALLED`). Эксперимент `aug34954` (эталон есть). Нужны QLKNN, TORBEAM, RABBIT.

### 2.4. Интегрированное моделирование «ядро + пьедестал»

**`imep_pw04`**, **`imep_pw08`** — «Template for ASTRA first approach to INTEGRATED MODELLING, CORE + EDGE» (Teobaldo Luda, 2018). Файлы отличаются только шириной пьедестала `CDWM1` = 0.04 / 0.08 (отсюда суффикс `pw04`/`pw08`). `IPEQL=5`, сетка 51×51. Решаются `NE`, `TE`, `TI` (`:EQ[2,CMHD4]`, `CMHD4=1` — до края) и ток `CU:EQ` (`MU=MUX`, `CU=CUX` как начальные). Транспорт — полуэмпирический: экспериментальные `HEXP`/`XEXP` в ядре, ступенчатые профили χ в пьедестале (`XSTEP`), адаптация ширины по `alfs` (`sbr/alfs.f90`), сглаживание `SMEARR`, плюс неоклассика NCLASS (`NEOCL4<::::;`, `sbr/nclass_mod.f90`, открытый) — она же даёт бутстреп `jbs_nc` и проводимость `cc_nc`. Пилообразные колебания — `mixint(CF2, CF6)::::;`. Источники и Zeff берутся из эксперимента: `ZEF=ZRD1X`, `NNCL=ZRD40X/1000`, `PEBM=car20x`, `PIBM=car21x`, `SNEBM=CAR18X`, `NIBM=car19x`, `PRAD=PRADX`. Набор `ZRD1`, `car17…car21`, `PRAD`, `MU`, `CU`, `BND` присутствует в `exp/30000_3.4` (DTT) — предположительно модели предназначены для него (не проверено; `ZRD40` там нет, тогда `NNCL=0`). Закрытых кодов не вызывают. Замечание: в блоке «Power block» есть строка `------------------------------------` без `!` в начале — как её обработает парсер, не проверено.

### 2.5. Стеллараторные и VMEC-модели

**`stelxes`** — тест стеллараторного обобщения уравнения тока на токамачном равновесии: `IPEQL=4` (SPIDER), внешняя индуктивность `lexthrs(ZRD96)<::::;` (Hirshman–Neilson), `CU:EQ` с `CC=CNSA`, `CD=0`, `CU=CC` в начале; профили заданы, `ZEF=2`. На графиках метрика (`G11`, `G22`, `G33`, `SG11`, `SG12`, `SG21`, `SG22`, `MV`) и баланс `abs(SG11*MV+SG12)`. `equ/log/stelxes` — в старом формате («Start file for version 7.1.0»). Эксперимент не установлен. Нужен SPIDER.

**`transpstell`** — стелларатор с VMEC: `IPEQL=7`, `vmec_option=0`, сетка 99×51. Профили и ток заданы; все транспортные коэффициенты обнулены; уравнение `F1:EQ` для радиального поля `ER=F1*1e3` с источником `SF1` из `dkes2astra_variables(j,4)` (модуль `src/astra/a2dkes.f90`). Вызов `DKES_INTERFACE_WRAPPER` в этой модели закомментирован. Предположительно для эксперимента `steltest` (W7-X, `exp/nml/steltest`: `&vmec`, `&travis` с `vmec_io/W7X-nc.input`). Нужны STELLOPT/VMEC (`$STELLOPT_HOME`) и таблицы DKES.

**`transpstell2`** — развитие `transpstell`: `vmec_option=9`, `DKES_INTERFACE_WRAPPER>::2.*TAU::;` (`sbr/dkes_wrapper.f90`, магистерская работа E. Buglione-Ceresa, 2025) — из DKES берутся `HE`, `XI`, проводимость `CC` и бутстреп `CD`, `ER` из амбиполярности. Эксперимент не установлен (`steltest` / `steltest2`). Нужны VMEC и таблицы DKES (namelist `ASTRA_DKES_INTERFACE`).

**`tok2vmec`** — токамак через VMEC: `IPEQL=7`, `vmec_option=4` (без Boozer), аналитические профили (`TE=0.0001*(1-XRHO**2)+0.1` и т. п.), `CU:EQ`. Вызывает `TOKSTELOUT>::::;` — такая подпрограмма **не найдена** нигде в дереве (`src/`, `sbr/`, `pyparse/`), т. е. модель не слинкуется. Предположительно для `exp/nml/tokstel` (файла `exp/tokstel` нет).

### 2.6. Прочие модели

**`fbe`** — free-boundary evolution на FEQIS: `IPEQL=5`, `ITFBE=2.5` (до t=2.5 с граница задана, после — свободная граница; `src/astra/stepup.F90`), `ICIRCQ=1` (решаются уравнения цепей катушек; до `ITFBE` токи катушек `CCOIL`/`VCOIL` берутся из эксперимента), `IPCTRL=0`, `ATREQ=0.005`, шаг 2 мс. Профили заданы, решается `CU:EQ`. NBI, RABBIT, TORBEAM закомментированы. Эксперимент `AUG33040_2500` (токи 12 катушек на t=2.5 с; эталон есть). Работает на Ubuntu без правок; расхождение с эталоном ~0.1 % по `CU`/`MU`, последний шаг (t=2.62) обрывается с `Error in find new boundary`, как и в эталоне.

**`astrahl_simple_msp`** — «Simple use of STRAHL-ASTRA8 with new output structure» (D. Fajardo, 2023). `IPEQL=4` (SPIDER). Три примеси: B (`AIM1=10`), W (184), Ar (40); их плотности, средние заряды, `ZEF` и `PRAD` считает STRAHL (`A2STRAHL(0, 1e6, 1e6, CSOL1)>::::;`, `sbr/a2strahl.f90`). Коэффициенты переноса примесей заданы аналитически (`CAR40…CAR45`), источники на краю `CPEL1..3 = CIMP1..3*1e19` ат/с. Для сравнения те же примеси решаются встроенными уравнениями `F1…F3:EQ`. Основные профили и ток заданы. В шапке модели написано «run e.g. with exp/AUG33040_2500», но эталон — `aug34954astrahl_simple_msp.CDF`. Нужны SPIDER и бинарники STRAHL (`$ASTRA_EXT/strahl/sep23/bin/strahl`, `result_to_astra`); атомные данные лежат в `astra-src/strahl/`.

**`showdata`** — демонстрация способов отображения данных; входные данные задаёт `exp/readme`. `NEQUIL=0`, `TE:AS`, `CU:AS` с аналитическим `MU`, сравнение профилей вдоль разных хорд (`CAR2X…CAR4X`). Написана в формате ASTRA 7 (строки-пояснения без `!`). Не запускается: нет `equ/log/showdata` (фатальная ошибка в `read_input.f90`), `exp/readme` ссылается на отсутствующие `udb/00000.tes`, `00000.tis`, `1D.tmp`, `2D.tmp`, а после патча парсера U-строк падает `AssertionError` в `parse_exp2d`.

### 2.7. Файлы `equ/log/`

`equ/log/<m>` — стартовые значения переменных и констант, которые ядро читает при запуске (`assign_val` в `read_input.f90`) и которые можно менять в окне (клавиша «I», `src/graph/gui_interaction.f90`). Формат: блоки `Variables:` (37 переменных: `AB`, `ABC`, `AIM1…3`, `AMJ`, `BTOR`, `ELONG`, `IPL`, `NNCL`, `QNBI`, `RTOR`, `SHIFT`, `TRIAN`, …), `Constants:` (`CF1…`, `CIMP1…`, `CSOL1…`, `CDHJ1…` и т. д.) и `Control parameters` (`DPOUT`, `TIME`, `TAUMIN`, `TAUMAX`, `TAUINC`, `DELVAR`, `TSCALE`, `NB2EQL`, `DTEQL`, …). 23 файла по 214 строк, `tok2vmec`, `transpstell`, `transpstell2` — 175–181 строка, `stelxes` — 175 строк в формате «Start file for version 7.1.0».

Совпадающие по содержимому группы (md5):

| одинаковые `log`-файлы |
|---|
| `flux_cuas`, `flux_emeq`, `flux_feqis`, `flux_nbi`, `flux_spider`, `flux_toric` |
| `flux_neo`, `flux_neo_tglf`, `flux_tglf` |
| `imep_pw04`, `imep_pw08`, `tglf` |
| `qlk`, `qlknn` |
| `fbe`, `flux_pel` |

`fbe` отличается от `flux_feqis` только `DPOUT` (0.01 против 0.2). У `astrahl_simple_msp` в `log` заданы источники примесей `CIMP1=20` (B), `CIMP2=0.5` (W), `CIMP3=5` (Ar).

### 2.8. Сводная таблица моделей

«Ubuntu» — запускается ли модель после нашего установщика (без закрытых кодов IPP и без TGLF/NEO/QuaLiKiZ/QLKNN/SPIDER/STRAHL/VMEC). «да*» — после комментирования строк `TORBA`/`RABBIT`. Всё, что не помечено «проверено», — вывод из чтения кода.

| модель | `IPEQL` | решаемые уравнения | внешние модули | эксперимент | Ubuntu |
|---|---|---|---|---|---|
| `fbe` | 5 FEQIS, свободная граница | `CU` | — (NEUT) | `AUG33040_2500` | да (проверено) |
| `flux_feqis` | 5 FEQIS | `CU` | TORBEAM, RABBIT | `aug34954`, `aug34954_t` | да* (проверено) |
| `flux_cuas` | 5 | — (`CU:AS`) | TORBEAM, RABBIT, EQDSK | `aug34954` | вероятно да* |
| `flux_emeq` | 1 EMEQ | `CU` | TORBEAM, RABBIT | `aug34954` (по аналогии) | вероятно да* |
| `flux_nbi` | 5 | `CU` | NBI (откр.), TORBEAM, EQDSK | `aug34954` | вероятно да* |
| `flux_pel` | 5 | `CU` | ABLATION, TORBEAM, RABBIT | не установлен | нет (нет данных `&pellet`) |
| `flux_spider` | 4 SPIDER | `CU` | SPIDER, TORBEAM, RABBIT | `aug34954` | нет |
| `flux_neo` | 5 | `CU` | NEO (IPC) | `aug34954` (по аналогии) | нет |
| `flux_neo_tglf` | 4 SPIDER | `CU` | SPIDER, TGLF, NEO | `aug34954` | нет |
| `flux_tglf` | 5 | `CU` | TGLF (IPC) | `aug34954` (по аналогии) | нет |
| `flux_tglf_serial` | 5 | `CU` | TGLF (lib) | `aug34954` (по аналогии) | нет |
| `flux_toric` | 5 | `CU` | TORIC, TORBEAM, RABBIT | `aug34954` | нет |
| `tglf` | 5 | `NE`, `TE`, `TI`, `CU` | TGLF, TORBEAM, RABBIT | `aug34954` | нет |
| `tglf_pid` | 5 | `NE`, `TE`, `TI`, `CU` | TGLF, TORBEAM, RABBIT | `AUG36982_3400` | нет |
| `tglf_rot_pred_08` | 5 | `NE`, `TE`, `TI`, `CU`, `UPAR` | TGLF, TORBEAM, RABBIT | не установлен | нет |
| `tglf_serial` | 5 | `NE`, `TE`, `TI`, `CU` | TGLF (lib), TORBEAM, RABBIT | не установлен | нет |
| `impli` | 5 (умолч.) | `NE`, `TE`+`TI` неявно, `CU` | TGLF, TORBEAM, RABBIT | не установлен | нет |
| `qlk` | 5 | `NE`, `TE`, `TI` | QuaLiKiZ, TORBEAM, RABBIT | `aug34954` | нет |
| `qlknn` | 5 | `NE`, `TE`, `TI` | QLKNN, TORBEAM, RABBIT | `aug34954` | нет |
| `imep_pw04` | 5 | `NE`, `TE`, `TI`, `CU` | NCLASS (откр.) | предположительно `30000_3.4` | неизвестно |
| `imep_pw08` | 5 | `NE`, `TE`, `TI`, `CU` | NCLASS (откр.) | предположительно `30000_3.4` | неизвестно |
| `astrahl_simple_msp` | 4 SPIDER | `F1…F3` (примеси) | STRAHL, SPIDER | `aug34954` / `AUG33040_2500` | нет |
| `stelxes` | 4 SPIDER | `CU` | SPIDER | не установлен | нет |
| `transpstell` | 7 VMEC | `F1` (`ER`) | VMEC/STELLOPT, DKES-данные | вероятно `steltest` | нет |
| `transpstell2` | 7 VMEC | `F1` (`ER`) | VMEC/STELLOPT, DKES | `steltest`/`steltest2` | нет |
| `tok2vmec` | 7 VMEC | `CU` | VMEC, `TOKSTELOUT` (не существует) | `tokstel` (только nml) | нет |
| `showdata` | — (`NEQUIL=0`) | — | — | `readme` | нет (проверено) |

Итого: на чистой Ubuntu реально доступны `fbe`, `flux_feqis` (проверено) и, вероятно, `flux_cuas`, `flux_emeq`, `flux_nbi` (после комментирования `TORBA`/`RABBIT`); `imep_pw04/08` закрытых кодов не требуют, но не проверены.

Важно для запуска на другой машине: FEQIS берёт геометрию катушек и стенки из `exp/cnf/<machine>_description_in.json`, где `<machine>` — ключ `-dev` у `as_exe` (по умолчанию `aug`). Для «Глобус-М2» понадобится свой JSON-файл (см. раздел 3.2).

---

## 3. Эксперименты `exp/`

### 3.1. Файлы экспериментов

| файл | машина / разряд | содержание |
|---|---|---|
| `aug34954` | ASDEX Upgrade #34954 | `NA1=91`, `RTOR=1.65`, геометрия камеры (`AWALL`, `ELONM`, `TRICH`, `AB`), D (`AMJ=2`). `IPL`, `BTOR` — U-файлы `MAG34954.*_AVG`; `TI`, `TE`, `NE`, `VTOR`, `CAR2` (q) — U-файлы из `udb/34954/` с множителями (`1e-3` эВ→кэВ, `1e-13` см⁻³→10¹⁹ м⁻³, `1e-2` см/с→м/с). Токи 12 катушек `CCOIL` на t=0, граница `BNDU` из `34954/34954.bnd.txt` (`_r`/`_z`). `ZRD3=ZRD8=2.5` (назначение не установлено). `FILTER 0.002`. Используется в 9 эталонах. |
| `aug34954_t` | ASDEX Upgrade #34954 | То же, но `TE` из нестационарного U-файла `E34954.IDA` (два среза: 4.0 и 4.8 с) вместо усреднённого, без `ZRD3`/`ZRD8`. Эталон `aug34954_tflux_feqis.CDF`. |
| `AUG33040_2500` | ASDEX Upgrade #33040, t≈2.5 с | `IPL=1.0`, `BTOR=2.5` заданы числами (U-файлы `MAG33040.*` закомментированы). Токи (`CCOIL`) и напряжения (`VCOIL`, нули) 12 катушек на t=2.5 с — нужны для свободной границы в `fbe`. Профили `TI`, `TE`, `NE`, `VTOR` (из `ANGF…`, т. е. угловая скорость), `CAR2` (q) и граница — из `udb/33040_2500/`. Эталон `AUG33040_2500fbe.CDF`. |
| `AUG36982_3400` | ASDEX Upgrade #36982, t=3.4 с | `NA1=161`; `RTOR`, `BTOR`, `IPL` на t=3.4. Профили `TEX`, `TIX`, `NEX`, `VTORX` (101 точка, `GRIDTYPE 12`), `CAR2` (80 точек) и граница `BND` (221 точка) — записаны прямо в файле («from DB: DB diag, fitted», усреднение 3.25–3.55 с). Эталон `AUG36982_3400tglf_pid.CDF`. |
| `30000_3.4` | DTT («dtt #30000 t=3.4», Divertor Tokamak Test) | `NA1=601`, `RTOR=2.19`, `AB=0.75`, `IPL=5.5`, `BTOR=5.85`, `ZRD1=1`, `ZRD2=46`. Профили `NE`, `TE`, `TI` (201 точка), `car17…car21`, `PRAD`, `MU`, `CU`, граница `BND` (92 точки). По набору величин подходит к `imep_pw04/08` (предположение). Своего `exp/nml/30000_3.4` нет. |
| `steltest` | шапка «aug #34954», но `RTOR=5.5` (масштаб стелларатора) | Профили `NEX`, `TEX` (91 точка, строки по 2186 символов) прямо в файле, `VTOR` и `CAR2` из U-файлов #34954, токи катушек AUG, `BNDU` закомментирован. Используется со стеллараторными моделями (`nml/steltest` содержит `&vmec`, `&travis` для W7-X). |
| `readme` | «ASDEX-Upgrade», учебный | Пример и описание формата файла эксперимента ASTRA 7: старый и новый ввод массивов (`POINTS`, `GRIDTYPE`, `NAMEXP`, `NTIMES`), разные типы сеток (0, 1, 2, 10, 12, 18, 19, 20 — хорды), ввод из U-файлов; после `END` — текст «Main Rules for Array Data Input». Ссылается на отсутствующие `00000.tes`, `00000.tis`, `1D.tmp`, `2D.tmp`. Пара к `equ/showdata`. Не запускается в A8 (см. 2.6). |

Подробно о формате — [04-experiment-data.md](04-experiment-data.md).

### 3.2. `exp/cnf/` — описание машины для FEQIS

`src/astra/machine_config.f90` читает `exp/cnf/<machine>_description_in.json`; `exe/greenMatrices.py` (вызывается из `as_exe.py` при каждом запуске) для каждого `*_description_in.json` пересчитывает `*_description_full.json` с матрицами Грина, если вход новее выхода. `src/feqis/feqis_solvers.f90` читает `_description_full.json`.

| файл | содержимое |
|---|---|
| `aug_description_in.json` | ASDEX Upgrade, полное описание: прямоугольная сетка FEQIS (`nR=nZ=65`, `Rmin=1.0`, `Rmax=2.3`, `Zmin=-1.6`, `Zmax=1.4`), `alpsep`, параметры бланкета (`res_blan`, `width_blan`), лимитер (`lim_*`, `Rlim`/`Zlim`, 113 точек), 21 элемент катушек (`R_coil`, `Z_coil`, `dR_coil`, `dZ_coil`, углы, число витков `m_turns`, группировка `m_equiv` в 12 цепей, `n_elem_coil`), матрица сопротивлений цепей 12×12 `resConduc`, пассивные элементы `x1…x4`, контуры камеры (`Rvessel`/`Zvessel`, 487 точек, 30 контуров). |
| `iter_description_in.json` | ИТЭР: только контуры для рисования (15 контуров, 397 точек); катушек нет — матрицы Грина не считаются (`greenMatrices.py` пишет «"Rmin" not found, skipping»). |
| `jet_description_in.json` | JET: один контур стенки (251 точка), катушек нет. |
| `read_json.py` | Скрипт matplotlib: рисует стенку и катушки из JSON (по умолчанию `aug_description_in.json`). |

Для расчёта со свободной границей на «Глобус-М2» нужен файл `exp/cnf/<имя>_description_in.json` в формате AUG и запуск с `-dev <имя>`.

### 3.3. `exp/nbi/`

`nbi/aug34954` — входной файл пакета NBI (Полевой) для 8 инжекторов AUG: для каждого источника строка `ZRDn` и 3 строки параметров (энергия, доли компонент, геометрия). Ищется ядром как `NBFILE` (`inquire_fname('nbi', exp_file, …)` в `read_input.f90`). Нужен модели `flux_nbi`.

### 3.4. `exp/nml/` — namelist-файлы внешних модулей

Файл `exp/nml/<эксперимент>` читается модулями по своим группам. Группы:

| группа | кто читает | смысл |
|---|---|---|
| `&rabbit_beam_geo`, `&partmix`, `&nbi_par`, `&species`, `&output_settings`, `&physics`, `&numerics` | RABBIT (`sbr/rabbit.F90`) | геометрия пучков, доли энергий, мощность (`pinj_file`), примесь, лимитер (`limiter_file`), таблицы |
| `&torbeam` | TORBEAM (`sbr/torba.F90`) | мощность и углы ECRH (`pecr_file`, `theta_file`, `phi_file`) |
| `&spider` | SPIDER (`src/astra/gs_solver.F90`) | `kprs`, `k_grids`, `epsros`, `enelss`, `key_plcs` |
| `&equil_settings` | FEQIS / `gs_solver` | допуски, релаксация, метод интерполяции, `max_iter` и т. п. |
| `&strahl_par` | STRAHL (`sbr/a2strahl.f90`) | шаг `tau_strahl`, координата, спад `ne`/`te` в SOL |
| `&toric` | TORIC (`sbr/torfpql_mod.F90`) | частота `freqcy=36.5e6`, цель по мощности, `nphi` |
| `&travis` | TRAVIS (`sbr/a2travis.f90`) | входной файл ECRH для стелларатора, номера пучков |
| `&vmec` | VMEC (`src/astra/a2vmec.f90`) | `mboz`, `nboz`, `phi_full_surfaces`, `dt_boozer`, `vac_phase_stel`, `f_vmec_settings` |
| `&pellet_ngs` | `sbr/ablation_ngs.f90` | модель абляции пеллеты для стелларатора (NGS) |

| файл | эксперимент | группы | ссылки на данные |
|---|---|---|---|
| `aug34954` | `aug34954` | RABBIT, `&torbeam`, `&spider`, `&equil_settings`, `&strahl_par`, `&toric` | `udb/34954/P34954.NBI_AVG`, `P34954.ECH_AVG`, `THE…`, `PHI…`, `rabbit/limiters/limiter_aug.dat` |
| `aug34954_t` | `aug34954_t` | **побайтно совпадает** с `aug34954` | те же |
| `AUG33040_2500` | `AUG33040_2500` | как `aug34954`, без `&toric` | `udb/33040_2500/*_AVG0` |
| `AUG36982_3400` | `AUG36982_3400` | как `aug34954`, без `&toric` | `udb/36982_3400/*` |
| `steltest` | `steltest` | RABBIT, `&torbeam`, `&travis` (W7-X), `&vmec`, `&strahl_par`, `&equil_settings`, `&spider`, `&pellet_ngs` | `udb/33040_2500/*`, `vmec_io/W7X-nc.input`, `vmec/templates/stell_files.nml` |
| `steltest2` | нет файла `exp/steltest2` | как `steltest` без `&pellet_ngs`; `&travis` с `vmec_io/WAL-nc.input` | `rabbit/limiters/limiter_walpha.dat` (в `astra-src/rabbit/limiters/` отсутствует) |
| `tokstel` | нет файла `exp/tokstel` | почти как `steltest2` (отличие — `&vmec`: `vac_phase_stel=1`, `phi_full_surfaces=0`) | то же |
| `dtt_exp_R2d19` | нет файла `exp/dtt_exp_R2d19` | RABBIT (2 пучка DTT, 510 кэВ), `&torbeam` | `udb/dtt/NBIdtt.W`, `udb/dtt/P34913.ECH_AVG`, `udb/dtt_test/ECRH30000.*` — в репозитории отсутствуют; `limiter_dtt.dat` есть |

### 3.5. `exp/icrh/`

`icrh/torica.inp` — шаблон входного файла TORIC 7.2 (namelist-ы `&toric_mode` с `toricmode='toric_fpql'`, `&loopdata` и др.). `sbr/torfpql_mod.F90` читает его из `exp/icrh/torica.inp`, подставляет аргументы вызова и пишет рабочий `torica.inp`. Нужен модели `flux_toric`.

---

## 4. Данные `udb/`

Все файлы — ASCII U-файлы (формат TRANSP/UFILES, ASDEX Upgrade, `AUGD`). «Ось Y» — радиальная координата или номер канала. Подробно о формате — [04-experiment-data.md](04-experiment-data.md).

### 4.1. `udb/34954/` (AUG #34954, для `aug34954`, `aug34954_t`, `steltest`)

| файл | тип | величина, единицы | срезы × точки | источник (из комментариев) |
|---|---|---|---|---|
| `34954.bnd.txt_r` | 2D | R границы, м | 10 × 180 | `write_u.py` |
| `34954.bnd.txt_z` | 2D | Z границы (в заголовке ошибочно «R BOUNDARY») | 10 × 180 | `write_u.py` |
| `E34954.IDA` | 2D | Te, эВ, по `rho_tor` | 2 (4.0; 4.8 с) × 51 | IDA, отображение EQH |
| `E34954.IDA_AVG` | 2D | Te, эВ | 1 × 51 | IDA, усреднение 2.00–2.94 с |
| `N34954.IDA_AVG` | 2D | ne, см⁻³ | 1 × 51 | IDA, 2.00–2.94 с |
| `I34954.CEZCMZ_AVG` | 2D | Ti, эВ | 1 × 51 | CEZ+CMZ, 2.00–2.85 с |
| `V34954.CEZCMZ_AVG` | 2D | Vφ, см/с | 1 × 51 | CMZ, 2.00–2.93 с |
| `Q34954.QPR` | 2D | q | 1 × 40 | EQI |
| `MAG34954.IPL_AVG` | 1D | ток плазмы, подпись «Amps», значение 0.797 (фактически МА) | 1 | `write_u.py` |
| `MAG34954.BTOR_AVG` | 1D | поле, подпись «Rp*Bt T.cm», значение 2.596 (фактически Тл) | 1 | `write_u.py` |
| `P34954.NBI_AVG` | 2D | мощность NBI по 8 источникам, Вт (работают №3 и №8, ~2.4 и 2.5 МВт) | 1 × 8 | `write_u.py` |
| `P34954.ECH_AVG` | 2D | мощность ECRH по 8 гиротронам, Вт (работает №1, 0.66 МВт) | 1 × 8 | `write_u.py` |
| `THE34954.ECH_AVG` | 2D | полоидальный угол ECRH, град | 1 × 8 | `write_u.py` |
| `PHI34954.ECH_AVG` | 2D | тороидальный угол ECRH, град | 1 × 8 | `write_u.py` |

### 4.2. `udb/33040_2500/` (AUG #33040, для `AUG33040_2500`, `steltest`)

Профили усреднены по 2.5006–3.4506 с.

| файл | тип | величина | срезы × точки | источник |
|---|---|---|---|---|
| `33040.bnd.txt_r`, `33040.bnd.txt_z` | 2D | R, Z границы, м | 20 × 90 | — |
| `TE33040.IDA_AVG0` | 2D | Te, эВ | 1 × 51 | IDA, EQH |
| `NE33040.IDA_AVG0` | 2D | ne, см⁻³ | 1 × 51 | IDA, EQH |
| `TI33040.CEZCMZ_AVG0` | 2D | Ti, эВ | 1 × 51 | CEZ+CMZ |
| `ANGF33040.CEZCMZ_AVG0` | 2D | Ω, рад/с | 1 × 51 | CEZ+CMZ `vrot` |
| `Q33040.QPR` | 2D | q | 1 × 40 | EQH |
| `MAG33040.IPL_AVG0` | 1D | ток, МА (0.78) | 1 | FPC `IpiFP`, 0–8 с (в `exp` закомментирован) |
| `MAG33040.BTOR_AVG0` | 1D | Btor, Тл (2.48) | 1 | FPC `BTF` (в `exp` закомментирован) |
| `Z33040.ZEF_AVG0` | 1D | Zeff (1.24) | 1 | не используется ни одним `exp` |
| `P33040.NBI_AVG0` | 2D | мощность NBI по 8 источникам, Вт | 1 × 8 | — |
| `P33040.ECH_AVG0` | 2D | мощность ECRH по 8 гиротронам, Вт (№5 0.63 МВт, №8 0.28 МВт) | 1 × 8 | — |
| `THEAS33040.ECH_AVG0`, `PHIAS33040.ECH_AVG0` | 2D | полоидальный / тороидальный угол ECRH, град | 1 × 8 | — |

### 4.3. `udb/36982_3400/` (AUG #36982, для `AUG36982_3400`)

Только данные нагрева (профили записаны прямо в `exp/AUG36982_3400`); все — `write_u.py`, t≈3.4 с.

| файл | величина |
|---|---|
| `P36982_3400.NBI_AVG` | мощность NBI по 8 источникам, Вт (ненулевой только №3, значение 10.06 — порядок величины в файле подозрительно мал для Вт) |
| `P36982_3400.ECH_AVG` | мощность ECRH, Вт (№1 0.60 МВт, №2 0.58 МВт) |
| `THE36982_3400.ECH_AVG`, `PHI36982_3400.ECH_AVG` | полоидальный и тороидальный углы ECRH, град |

Файлы, на которые ссылаются `exp/nml/dtt_exp_R2d19` (`udb/dtt/`, `udb/dtt_test/`) и `exp/readme` (`udb/00000.tes`, `00000.tis`), в репозитории отсутствуют.

---

## 5. Формулы `fml/`

Формула — фрагмент Fortran, который парсер вставляет в генерируемый код внутри цикла по радиусу `J` (см. [03-model-language.md](03-model-language.md)); имя файла в нижнем регистре, в модели — в любом регистре. Формулы могут включать друг друга (`include 'fml/…'`). Колонка «§4»: «да» — величина перечислена в `docu/section4.tex`, «List of expressions».

Общие замечания по каталогу:

- 12 формул вызывают подпрограммы, которых в A8 нет нигде в дереве: `MIX`, `MIX1` (`ccmhd`, `ccmhd1`, `hmhd1`, `hmhd2`, `xmhd1`, `xmhd2`) и `SMOFML` (`habms`, `haets`, `hagbs`, `harls`, `hatls`, `hetis`). Модель с такой формулой, скорее всего, не слинкуется (не проверено сборкой).
- `eli3` вызывает функцию `elin3R`, которая в дереве не найдена.
- `pdtf` в шапке помечена как устаревшая (Pereverzev, 2002): «cannot be called from a model any more».
- В нескольких файлах заголовок скопирован из другого файла и не соответствует имени: `zikr` и `prkr` (заголовок «Argon», по имени — криптон), `zineo`, `zinit`, `ziwol` (заголовок «Average C charge», по имени — Ne, N, W), `zine` («Average Ne chNEge»), `peht`/`piht` (заголовки `PEDT`/`PIDT`), `peiad` (заголовок `PEIAU`), `hatly` (заголовок `HATLI`), `dchr2`/`hchr2`/`xchr2` (заголовки `DCHR`/`HCHR`/`XCHR`), `betpt` (заголовок `betpl`).

### 5.1. Радиальные функции и сетка

| файл | описание | §4 |
|---|---|---|
| `fa` | малый радиус в экваториальной плоскости (то же, что `AMETR`), м | да |
| `flin` | линейный профиль | да |
| `fpa` | параболический профиль по «a» | да |
| `fpr` | параболический профиль | да |
| `fr` | радиус, м | да |
| `frs` | радиус во вспомогательных узлах сетки, м | да |
| `fx` | нормированный радиус r/r_edge (то же, что `FLIN`) | да |
| `rnode` | номер радиального узла, ближайшего к `RHO(j)` | — |

### 5.2. Проводимость и эффективность генерации тока

| файл | описание | §4 |
|---|---|---|
| `ccneu` | классическая проводимость, МА/(В·м) | да |
| `ccsp` | спитцеровская проводимость | да |
| `ccspx` | спитцеровская проводимость (вариант) | да |
| `ccspsa` | спитцеровская проводимость для модели Sauter | — |
| `chotf` | проводимость «hot», Fisch | да |
| `cnhh` | неоклассическая проводимость Hinton–Hazeltine | да |
| `cnhr` | неоклассическая проводимость Hirshman | да |
| `cnsa` | неоклассическая проводимость Sauter–Angioni–Lin-Liu (PoP 6, 2834, 1999) | — |
| `ccmhd` | аномальная проводимость из-за пилообразных колебаний (q<1); вызывает `MIX` (нет в A8) | да |
| `ccmhd1` | то же, с искусственным удержанием q>1 (2006); вызывает `MIX1` (нет в A8) | — |
| `cuohm` | плотность омического тока, МА/м² | да |
| `eflhn` | эффективность генерации тока НГ-волной, узкий спектр, м/В | да |
| `eflhw` | эффективность генерации тока НГ-волной, широкий спектр | да |
| `fte` | множитель без описания в заголовке; пример: `CC=CCHR*min(CNHR/CCSP, FTE)` | да |

### 5.3. Бутстреп-ток и неоклассические коэффициенты токового уравнения

Тройки `DC`/`HC`/`XC` — коэффициенты при градиентах плотности, `Te` и `Ti` в бутстреп-токе (используются как `DC=DCKIM; HC=HCKIM; XC=XCKIM`).

| файл | описание | §4 |
|---|---|---|
| `dcha`, `hcha`, `xcha` | бутстреп-ток (модель в заголовке не указана) | да |
| `dchh`, `hchh`, `xchh` | Hinton–Hazeltine | да |
| `dchr`, `hchr`, `xchr` | Hirshman (Phys. Fluids 31, 3150, 1988) | да |
| `dchr2`, `hchr2`, `xchr2` | Hirshman с поправкой на столкновительность (Sauter 1994, Pacher 1997) | — |
| `dckim`, `hckim`, `xckim` | Kim (Phys. Fluids B 3, 2050, 1991), все режимы столкновительности | да |
| `dckm1`, `hckm1`, `xckm1` | Kim, вариант | да |
| `dhkim` | вспомогательный include-файл для `HCKIM`/`XCKIM` | да |
| `dcsa`, `hcsa`, `xcsa` | Sauter | — |
| `enhh0`, `enhh2` | «electrical heat convection», Hinton–Hazeltine | да |
| `enhh1` | «electrical particle convection», Hinton–Hazeltine | да |
| `vras` | неоклассический пинч Уэра | — |
| `vrhh` | пинч Уэра, Hinton–Hazeltine | — |

### 5.4. Неоклассический перенос

| файл | описание | §4 |
|---|---|---|
| `hnasi` | ионная теплопроводность Angioni–Sauter | — |
| `hnchi` | ионная теплопроводность Chang–Hinton | да |
| `hchgp` | ионная теплопроводность Chang–Hinton (вариант) | да |
| `hchii` | улучшенная Chang–Hinton | да |
| `xch86` | Chang–Hinton (1986) | да |
| `hngse` | электронная, Galeev–Sagdeev | да |
| `hngsi` | ионная, Galeev–Sagdeev | да |
| `hngsb` | ионная, Galeev–Sagdeev, банановый режим | да |
| `hngsp` | ионная, Galeev–Sagdeev, плато | да |
| `hnpsi` | ионная, Пфирш–Шлютер | да |
| `fowc` | поправка на конечную ширину орбит к неоклассической ионной теплопроводности (Lin, Tang, Lee 1997) | да |
| `dnneo` | неоклассическая диффузия частиц (диссертация Angioni) | — |
| `vpsww` | неоклассическая полоидальная скорость вращения, м/с | — |

### 5.5. Аномальный перенос: эмпирические и полуэмпирические модели

| файл | описание | §4 |
|---|---|---|
| `cerl` | тепловой пинч Rebut–Lallia, м/с | да |
| `dbohm` | бомовская диффузия | да |
| `haalc` | Alcator | да |
| `habm` | Bohm/gyroBohm, Taroni (L-мода) | да |
| `habms` | сглаженная `HABM`; вызывает `SMOFML` (нет в A8) | да |
| `habom` | бомовская теплопроводность | да |
| `hibom` | бомовская теплопроводность (заголовок как у `habom`; вероятно, ионная) | — |
| `hacte` | бесстолкновительные запертые электроны | да |
| `haed` | электронный дрейф | да |
| `haeti`, `hetai` | η_i-мода | да |
| `haets`, `hetis` | η_i сглаженная; вызывают `SMOFML` (нет в A8) | да |
| `hagb` | gyroBohm | да |
| `hagbs` | gyroBohm сглаженная; вызывает `SMOFML` | да |
| `hagbs1` | gyroBohm сглаженная, вариант | да |
| `haitf` | коэффициент для теплопроводности по скейлингу ИТЭР | да |
| `haitr` | теплопроводность по скейлингу ИТЭР | да |
| `hamm` | Мережкин–Муховатов | да |
| `hanag` | Neo-Alcator–Goldston | да |
| `hanal` | Neo-Alcator | да |
| `hapa` | Парайль | да |
| `hapue` | Парайль–Юшманов–Есипчук | да |
| `hapyu` | Парайль–Юшманов | да |
| `harl` | Rebut–Lallia | да |
| `harls` | Rebut–Lallia сглаженная; вызывает `SMOFML` | да |
| `harnq` | зависимость от ρ, ν, q | да |
| `harpl` | rippling-мода | да |
| `hascl` | скиновая, столкновительная | да |
| `hasdl` | по скейлингу L-моды ASDEX (Karulin 1993) | — |
| `hatl` | Taroni (L-мода) | да |
| `hatli` | Taroni (L-мода), вариант | да |
| `hatly` | ещё один вариант Taroni (заголовок `HATLI`) | — |
| `hatls` | Taroni сглаженная; вызывает `SMOFML` | да |
| `hbjet` | бомовский член транспортной модели JET | да |
| `hgbej` | gyroBohm-член модели JET (электроны) | да |
| `hgbij` | gyroBohm-член модели JET (ионы) | да |
| `hit89` | формула Юшманова для χ (ИТЭР-89) | да |
| `xirl` | ионная Rebut–Lallia, сглаженная | да |
| `tecrj` | критический градиент Te (Janeschitz) | — |
| `tecrl` | критический градиент Te (Rebut–Lallia) | да |
| `tecrl1` | то же, вариант | — |

### 5.6. Пилообразные колебания и МГД

| файл | описание | §4 |
|---|---|---|
| `almhd` | α_MHD(r) | да |
| `betpv` | полоидальная бета для модели пилы | — |
| `haq1`, `haq1c` | аномальная теплопроводность при q<1 | да |
| `hmhd1`, `hmhd2` | то же для `HE`; вызывают `MIX`/`MIX1` (нет в A8) | да |
| `xmhd1`, `xmhd2` | то же для `XI`; вызывают `MIX`/`MIX1` (нет в A8) | да |

### 5.7. Экспериментальные и эффективные коэффициенты (баланс мощности)

| файл | описание | §4 |
|---|---|---|
| `dndif`, `dnexp` | экспериментальный коэффициент диффузии частиц | да |
| `hexp` | экспериментальная электронная теплопроводность | да |
| `hegn` | экспериментальная электронная теплопроводность (вариант с конвекцией) | да |
| `heeff` | эффективная электронная теплопроводность | да |
| `xexp` | экспериментальная ионная теплопроводность | да |
| `xign` | экспериментальная ионная теплопроводность (вариант) | да |
| `xieff` | эффективная ионная теплопроводность | да |
| `xtexp` | эффективная диффузия тороидального импульса | — |

### 5.8. Микронеустойчивости, критические градиенты, сдвиг

| файл | описание | §4 |
|---|---|---|
| `gaitg` | линейный инкремент ITG | — |
| `gitg` | линейный инкремент ITG (из процедуры Dorland) | — |
| `gitg0` | множитель перед `(RTOR/LTI-RLTcr)` в инкременте ITG | — |
| `rltcr` | R/L_T,crit для ITG (модель IFS/PPPL) | да |
| `rltcz` | R/L_T,crit, углеродная ветвь (IFS/PPPL) | да |
| `rltkd` | R/L_T,crit для ITG (IFS/PPPL, вариант) | — |
| `rltwn` | R/L_T,crit для ITG (Weiland–Nordman) | да |
| `rlewn` | R/L_Te,crit для TEM (Weiland–Nordman) | — |
| `shat`, `sheas` | шир для модели IFS/PPPL | — |
| `d2ti` | ~сдвиг полоидальной скорости по Парайлю (PPCF 40, 805, 1998) | да |
| `rotsh` | скорость сдвига вращения | да |
| `rhos` | критическое ρ* (Emi, 2004) | — |

### 5.9. Вращение и скорости

| файл | описание | §4 |
|---|---|---|
| `cvfi` | V_φ = CVFI·χ_φ (пинч импульса), 1/м | — |
| `cxfi` | отношение χ_φ/χ_Ti (Angioni) | — |
| `vdia` | диамагнитная скорость | да |
| `vdie` | электронная диамагнитная скорость | — |
| `cs` | скорость звука √(Te/Mi) | да |
| `vsi` | ионно-звуковая скорость √((Te+Ti)/Mi) | да |
| `vte` | тепловая скорость электронов | да |
| `vti` | тепловая скорость ионов | да |

### 5.10. Термоядерные реакции и альфа-частицы

| файл | описание | §4 |
|---|---|---|
| `pdt` | мощность D-T синтеза, МВт/м³ | да |
| `pdtf` | D-T с учётом NBI; **устарела**, несовместима с текущей NBI | — |
| `pedt`, `pedt1` | доля мощности D-T синтеза на электроны | да |
| `peht` | то же (заголовок `PEDT`, через `paion1`, C. M. Roach 1999) | — |
| `pedtf` | на электроны с учётом быстрых ионов NBI | — |
| `pidt`, `pidt1` | доля мощности D-T синтеза на ионы | да |
| `piht` | то же (заголовок `PIDT`, через `paion1`) | — |
| `pidtf` | на ионы с учётом быстрых ионов | — |
| `paion`, `paion1` | доля мощности альфа-частиц, отдаваемая ионам | да / — |
| `qdtf` | полная мощность синтеза с быстрыми ионами, МВт | — |
| `qedt`, `qidt` | мощность синтеза на электроны / ионы, МВт | да |
| `qedtf`, `qidtf` | то же с быстрыми ионами | — |
| `nalph` | плотность быстрых альфа-частиц | — |
| `nuas` | частота торможения альфа-частиц D-T | — |
| `ts_alf` | время торможения альфа-частиц (Stix 1972) | — |
| `vca` | критическая скорость альфа-частиц | — |
| `presal` | давление альфа-частиц | — |
| `calpha` | аномальная конвекция альфа-частиц | — |
| `dalpha` | аномальная диффузия альфа-частиц | — |
| `cphe` | чисто конвективный коэффициент потока He | — |
| `cthe` | коэффициент термодиффузии He | — |
| `svdt` | скорость реакции D-T (аппроксимация) | да |
| `svd1`, `svd2`, `svdbh` | скорость реакции D-D | да |
| `svdhe` | скорость реакции D-³He | да |

### 5.11. Быстрые ионы NBI

| файл | описание | §4 |
|---|---|---|
| `ibm` | ток пучка внутри поверхности ρ, МА | да |
| `nufe` | частота торможения быстрых ионов | — |
| `pbeie` | потери электронов на ионизацию нейтралов пучка | да |
| `pbicx` | нагрев ионов перезарядкой на нейтралах пучка | да |
| `bpfit` | локальная полоидальная бета с быстрыми ионами | — |

### 5.12. Обмен энергией, омический нагрев, конвективные потоки

| файл | описание | §4 |
|---|---|---|
| `pei`, `pei2` | электрон-ионный обмен (коэффициент при ΔT), МВт/м³/кэВ | да / — |
| `peiau`, `peiad` | электрон-ионный обмен, варианты (`peiad` — заголовок `PEIAU`) | — |
| `peicl` | электрон-ионный обмен, МВт/м³ | да |
| `peicl2` | то же, вариант | — |
| `pehcl` | обмен электронов с ионами H | да |
| `phicl` | термализация H–ионы | да |
| `suzpei` | Σ Nz·Z²/Mz тепловых ионов для обмена энергией | — |
| `poh` | омический нагрев | да |
| `pohe` | омический нагрев «как должно быть» (фрикционный) | — |
| `pjoul` | джоулев нагрев | да |
| `pegn`, `pign` | мощность, переносимая потоком частиц (электроны / ионы) | да |
| `peign` | обмен энергией из-за потока частиц | да |
| `qeign` | ∫ `PEIGN` dV | да |
| `qign` | конвективный поток тепла ионов с частицами | да |

### 5.13. Нейтралы и атомные процессы

| файл | описание | §4 |
|---|---|---|
| `peneu` | потери электронов на ионизацию холодных нейтралов | да |
| `penli` | потери электронов на излучение нейтралов | да |
| `picx` | нагрев ионов перезарядкой на холодных нейтралах | да |
| `pineu` | источник тепла ионов от холодных нейтралов | да |
| `pionz` | источник тепла ионов от ионизации | да |
| `pirec` | потери на рекомбинацию | да |
| `pitcx` | источник тепла ионов от перезарядки | да |
| `snneu` | S/Ne от нейтралов (ионизация + перезарядка − рекомбинация) | да |
| `snnie` | S = ⟨σv⟩·Nn | да |
| `snnii` | S = ⟨σv⟩·Nn·Ni | да |
| `snnr` | S = −⟨σv⟩·Ne² (рекомбинация) | да |
| `svcx`, `svcxx` | скорость перезарядки (аппроксимации) | да |
| `svie`, `sviep`, `sviex` | ионизация электронным ударом | да |
| `svii` | ионизация ионным ударом | да |
| `svrec`, `svrecx` | рекомбинация | да |

### 5.14. Излучение

| файл | описание | §4 |
|---|---|---|
| `pbrad` | электронное тормозное излучение | да |
| `pbr1`, `pbr2`, `pbr3` | тормозное излучение 1-й/2-й/3-й примеси | — |
| `psync` | синхротронное излучение | да |
| `qbrem` | полная мощность тормозного излучения, МВт | да |
| `prarg` | функция излучения аргона, МВт·м³ | — |
| `prber` | бериллий (Z=4, A=9) | — |
| `prfer` | железо | да |
| `prkr` | заголовок «Argon» (Leonov 1999), по имени — криптон | — |
| `prneo` | неон | да |
| `prnit` | азот | да |
| `proxi`, `proxy` | кислород (два варианта) | да |
| `prwol` | вольфрам | да |
| `prwol_puet2010` | функция охлаждения вольфрама (Pütterich 2010) | — |
| `prxen` | ксенон | — |

### 5.15. Средний заряд примесей (корональное приближение)

| файл | описание | §4 |
|---|---|---|
| `zicar` | углерод | да |
| `ziar` | аргон | — |
| `zikr` | заголовок «Ar» (Pacher 1993), по имени — криптон | — |
| `zibe` | бериллий | — |
| `zifer` | железо | да |
| `zine` | неон | — |
| `zineo` | заголовок «C», по имени — неон | да |
| `zinit` | заголовок «C», по имени — азот | да |
| `zioxi` | кислород | да |
| `ziwol` | заголовок «C», по имени — вольфрам | да |
| `ziwo2` | средний квадрат заряда W | — |
| `ziwol_puet2010` | средний заряд W (Pütterich 2010) | — |
| `zixen` | ксенон | — |
| `zixe2` | средний квадрат заряда Xe | — |
| `zzef` | расчёт Zeff «according to guidelines» | — |

### 5.16. SOL

| файл | описание | §4 |
|---|---|---|
| `petsl` | продольные потери тепла электронов в SOL (Becker, NF 1995) | да |
| `pitsl` | то же для ионов | да |

### 5.17. Столкновения и частоты

| файл | описание | §4 |
|---|---|---|
| `coulg` | кулоновский логарифм | да |
| `clgesa`, `clgisa` | кулоновский логарифм для электронов / ионов (Sauter) | — |
| `nue` | полная частота электронных столкновений | да |
| `nuee` | частота e-e столкновений | да |
| `nues` | относительная частота электронных столкновений ν*_e | да |
| `nuessa` | ν*_e (Sauter) | — |
| `nui` | полная частота ионных столкновений | да |
| `nuis` | ν*_i | да |
| `nuissa` | ν*_i (Sauter) | — |
| `nupp` | частота p-p столкновений | да |
| `omde`, `omdi` | электронная / ионная диамагнитная частота | — |

### 5.18. Давление и бета

| файл | описание | §4 |
|---|---|---|
| `bete`, `beti` | локальная тороидальная бета электронов / ионов | да |
| `betpl` | локальная полоидальная бета (электроны) | да |
| `betpt` | локальная полоидальная бета полная (Fable 2014) | — |
| `prese`, `presi` | давление электронов / ионов, МДж/м³ | да |
| `prest` | полное давление | да |

### 5.19. Характерные длины и безразмерные градиенты

| файл | описание | §4 |
|---|---|---|
| `lne`, `lni`, `lnz1` | характерная длина профиля ne, ni, NIZ1 | да |
| `lte`, `lti` | характерная длина Te, Ti | да |
| `etae`, `etai` | d ln Te / d ln ne, d ln Ti / d ln ne | да |
| `etan` | (d ln ne)⁻¹/R | да |
| `rlne` | R/L_ne | — |
| `rlf1`, `rlf2` | R/L для `F1`, `F2` | — |
| `rli` | ларморовский радиус иона | да |
| `rls` | ρ* (по Te) | да |
| `sqz` | squeezing-фактор основных ионов | да |

### 5.20. Ток, q, индуктивность, поля

| файл | описание | §4 |
|---|---|---|
| `didt` | dI/dt полного тока, МА/с | да |
| `eli3` | li(3) из решателя равновесия; вызывает `elin3R` (не найдена) | — |
| `qbm`, `qbs`, `qcd`, `qnind` | q от тока NBI / бутстрепа / CD / неиндукционного тока | да |
| `edr` | поле Драйсера | да |
| `epar` | продольное электрическое поле | да |
| `epl` | тороидальное электрическое поле | да |

### 5.21. Запертые частицы

| файл | описание | §4 |
|---|---|---|
| `tef` | доля запертых электронов (банановый режим) | — |
| `tpf` | доля запертых частиц | да |

(`fte` — см. 5.2.)

### 5.22. Интегралы источников уравнений `F0…F9`, `NE`, `TE`, `TI`

| файл | описание | §4 |
|---|---|---|
| `qdv0`, `qdv1`, `qdv2`, `qdv3`, `qdv4`, `qdv5`, `qdv6`, `qdv7`, `qdv8`, `qdv9` | ∫₀ᴿ `SDV0…SDV9` dV (источники уравнений `F0…F9`; во всех заголовках единицы записаны как `[F0/s]`) | — |
| `qdve`, `qdvi` | ∫ `PDVE`, `PDVI` dV, МВт | — |
| `qdvn` | ∫ `SDVN` dV, 10¹⁹/с | — |

### 5.23. Скейлинги времени удержания

| файл | описание | §4 |
|---|---|---|
| `tau89` | ITER-89 (Юшманов) | да |
| `tauna` | Neo-Alcator для ИТЭР | да |
| `titer` | формула Юшманова для τE ИТЭР | да |
| `thq99` | τE,th IPB98(y,2), H-мода с ELM | — |
| `tlq97` | τE,th ITER L-мода (1997) | — |

---

## 6. Функции `fnc/`

Функции компилируются в `lib/libfnc.a` (`exe/Makesub sub=fnc`). Имя функции = имя выражения + `R`; аргумент — радиус (`YR`, м), до которого берётся интеграл/среднее, например `QETOT` → `QETOTR(YR)`; в модели обычно пишут `QETOTB` (на границе) или `QETOT(r)`. Все — `double precision function`. «§4» — упоминание в `docu/section4.tex`.

### 6.1. Интегралы мощности и частиц

| файл | сигнатура | описание | §4 |
|---|---|---|---|
| `qbrad.f90` | `QBRADR(YR)` | ∫ `PBRAD` dV, МВт | — |
| `qbtot.f90` | `QBTOTR(YR)` | ∫ `PBEAM` dV | да |
| `qdt.f90` | `QDTR(YR)` | ∫ `PDT` dV | да |
| `qde.f90`, `qdi.f90` | `QDER(YR)`, `QDIR(YR)` | ∫ `PDE`, `PDI` dV | — |
| `qdn.f90` | `QDNR(YR)` | ∫ `SDN` dV, 10¹⁹/с | — |
| `qd0.f90`, `qd1.f90`, `qd2.f90`, `qd3.f90`, `qd4.f90`, `qd5.f90`, `qd6.f90`, `qd7.f90`, `qd8.f90`, `qd9.f90` | `QD0R(YR)` … `QD9R(YR)` | ∫ `SD0…SD9` dV | — |
| `q0tot.f90`, `q1tot.f90`, `q2tot.f90`, `q3tot.f90`, `q4tot.f90`, `q5tot.f90`, `q6tot.f90`, `q7tot.f90`, `q8tot.f90`, `q9tot.f90` | `Q0TOTR(YR)` … `Q9TOTR(YR)` | ∫ `SF0TOT…SF9TOT` dV | — |
| `qebm.f90`, `qibm.f90` | `QEBMR(YR)`, `QIBMR(YR)` | ∫ `PEBM`, `PIBM` dV | — |
| `qedwt.f90`, `qidwt.f90` | `QEDWTR(YR)`, `QIDWTR(YR)` | dWe/dt, dWi/dt | да |
| `qegn.f90` | `QEGNR(YR)` | ∫ `PEGN` dV | да |
| `qeiad.f90`, `qeiau.f90` | `QEIADR(YR)`, `QEIAUR(YR)` | ∫ `PEIAD`, `PEIAU` dV | — |
| `qeicl.f90` | `QEICLR(YR)` | ∫ `PEICL` dV | да |
| `qeneu.f90`, `qineu.f90` | `QENEUR(YR)`, `QINEUR(YR)` | ∫ `PENEU`, `PINEU` dV | да |
| `qicx.f90` | `QICXR(YR)` | ∫ `PICX` dV | — |
| `qetot.f90`, `qitot.f90` | `QETOTR(YR)`, `QITOTR(YR)` | ∫ `PETOT`, `PITOT` dV | да |
| `qex.f90`, `qix.f90` | `QEXR(YR)`, `QIXR(YR)` | ∫ `PEX`, `PIX` dV | да |
| `qjoul.f90` | `QJOULR(YR)` | ∫ `PJOUL` dV | да |
| `qoh.f90` | `QOHR(YR)` | ∫ `POH` dV | да |
| `qrad.f90`, `qradx.f90` | `QRADR(YR)`, `QRADXR(YR)` | ∫ `PRAD`, `PRADX` dV | да |
| `qsync.f90` | `QSYNCR(YR)` | ∫ `PSYNC` dV | да |
| `qtot.f90` | `QTOTR(YR)` | ∫ (`PETOT`+`PITOT`) dV | да |
| `qtok.f90` | `QTOKR(YR)` | ∫ (`PETOT`+`PITOT`) dV (дубль `qtot` по описанию) | — |
| `qntot.f90` | `QNTOTR(YR)` | ∫ `SNTOT` dV, 10¹⁹/с | да |
| `qnx.f90` | `QNXR(YR)` | ∫ `SNX` dV | да |
| `qndnt.f90` | `QNDNTR(YR)` | d/dt ∫ `NE` dV | да |
| `qudut.f90` | `QUDUTR(YR)` | d/dt ∫ `Upar·Ups` dV | — |
| `qutot.f90` | `QUTOTR(YR)` | интеграл, связанный с моментом `TTRQ` (в заголовке только имя) | — |
| `shine.f90` | `SHINER(YR)` | ∫ `SVD1·NDEUT²` dV — выход DD-нейтронов (заголовок `SDDN`) | — |
| `ptot.f90` | `PTOTR(YR)` | полная мощность нагрева, МВт | — |

### 6.2. Средние, энергосодержание, число частиц

| файл | сигнатура | описание | §4 |
|---|---|---|---|
| `neav.f90`, `nexav.f90` | `NEAVR(YR)`, `NEXAVR(YR)` | средняя по объёму плотность (расчёт / эксперимент) | да |
| `nech.f90` | `NECHR(r_in)` | средняя по хорде плотность | да |
| `neca.f90` | `NECAR()` | простое среднее `NE` по узлам (заголовок «NECH») | — |
| `net.f90` | `NETR(YR)` | полное число электронов | — |
| `teav.f90`, `texav.f90` | `TEAVR(YR)`, `TEXAVR(YR)` | средняя по объёму Te (в заголовке ошибочно «ion») | да |
| `tiav.f90`, `tixav.f90` | `TIAVR(YR)`, `TIXAVR(YR)` | средняя по объёму Ti | да |
| `tendn.f90`, `tindn.f90` | `TENDNR(YR)`, `TINDNR(YR)` | средневзвешенная по плотности Te / Ti | да |
| `zndn.f90` | `ZNDNR(YR)` | средневзвешенный по плотности Zeff | да |
| `vol.f90` | `VOLR(YR)` | объём, м³ | — |
| `we.f90`, `wex.f90` | `WER(YR)`, `WEXR(YR)` | ∫ 3/2·ne·Te dV (расчёт / эксперимент), МДж | да |
| `wi.f90`, `wix.f90` | `WIR(YR)`, `WIXR(YR)` | ∫ 3/2·ni·Ti dV | да |
| `wtot.f90`, `wtotx.f90` | `WTOTR(YR)`, `WTOTXR(YR)` | полная тепловая энергия | да |
| `wtoz.f90` | `WTOZR(YR)` | энергия с быстрыми частицами (`PFAST`, `PBLON`, `PBPER`) | — |
| `wc.f90`, `wcx.f90` | `WCR(YR)`, `WCXR(YR)` | энергия ядра (без пьедестала), расчёт / эксперимент | — |
| `wce.f90`, `wcex.f90` | `WCER(YR)`, `WCEXR(YR)` | то же, электроны | — |
| `wci.f90`, `wcix.f90` | `WCIR(YR)`, `WCIXR(YR)` | то же, ионы | — |
| `wbpol.f90` | `WBPOLR(YR1)` | энергия полоидального поля, МДж | да |

### 6.3. Времена удержания

| файл | сигнатура | описание | §4 |
|---|---|---|---|
| `taue.f90` | `TAUER(YR)` | полное время удержания энергии | да |
| `tauee.f90`, `tauei.f90` | `TAUEER(YR)`, `TAUEIR(YR)` | электронное / ионное | да |
| `taup.f90` | `TAUPR(YR)` | время удержания частиц | да |
| `tauf0.f90`, `tauf1.f90`, `tauf2.f90`, `tauf3.f90`, `tauf4.f90`, `tauf5.f90`, `tauf6.f90`, `tauf7.f90`, `tauf8.f90`, `tauf9.f90` | `TAUF0R(YR)` … `TAUF9R(YR)` | время удержания величин `F0…F9` | — |
| `h98iter.f90` | `H98ITERR(YR)` | скейлинг IPB98(y,2) (комментарий: в τ_scal только P_abs, не P_sep) | — |

### 6.4. Ток, индуктивность, q, бета

| файл | сигнатура | описание | §4 |
|---|---|---|---|
| `ibs.f90` | `IBSR(YR)` | бутстреп-ток внутри r, МА | да |
| `icd.f90` | `ICDR(YR)` | генерируемый ток | да |
| `iohm.f90` | `IOHMR(YR)` | омический ток | да |
| `itot.f90` | `ITOTR(YR)` | полный ток (в заголовке ошибочно «IBS») | да |
| `li3.f90` | `LI3R(YRO)` | внутренняя индуктивность li(3) | — |
| `lint.f90` | `LINTR(YRO)` | внутренняя индуктивность li | да |
| `licd.f90` | `LICDR(YR1)` | li от доли генерируемого тока | да |
| `qmin.f90` | `QMINR(YR)` | min q на [0, a] | да |
| `xqmin.f90` | `XQMINR(YR)` | положение min q | да |
| `el095.f90` | `el095r(YD)` | вытянутость на поверхности 95 % потока; переменная `JJ1` не инициализирована перед циклом (по коду) | — |
| `ftllm.f90` | `FTLLMR(YR)` | эффективная доля запертых частиц (Lin-Liu, Miller 1995) | да |
| `beta.f90` | `BETAR(YR)` | тороидальная бета, % (с давлением NBI) | — |
| `betan.f90` | `BETANR(YR)` | нормированная бета βN (с NBI и быстрыми) | — |
| `betat.f90` | `BETATR(YR)` | нормированная тепловая βN | — |
| `betaj.f90` | `BETAJR(YR)` | полоидальная бета (r) | да |
| `betbm.f90` | `BETBMR(YR)` | полоидальная бета пучка | да |
| `betp3.f90` | `BETP3R(YR)` | βpol = Wp/Wm (2008) | — |
| `betr.f90` | `BETRR(YR)` | фактор Тройона β[%]/(I/aB) | да |
| `alim.f90` | `ALIMR(YR)` | α_MHD на пределе баллонной устойчивости | — |
| `prcar.f90` | `PRCARR(YR)` | функция излучения углерода, МВт·м³ | да |

---

## 7. Подпрограммы `sbr/`

Вызов из модели: `ИМЯ(аргументы)` с модификаторами времени (`<`, `>`, `:t:t0:t1:`), см. [03-model-language.md](03-model-language.md). `docu/section4.tex` («Plug-in subroutines») описывает `GNEX`, `GNXSRC`, `MINMAX`, `MIXINT`, `NEUT`, `SETNAV`, `SMEARR`, `STARR`/`STVAR`, `TSCTRL` и `NBI`; из них в A8 есть только `GNEX`, `MIXINT`, `NEUT`, `SMEARR`, `NBI` — остальных (`GNXSRC`, `MINMAX`, `SETNAV`, `STARR`, `STVAR`, `TSCTRL`) в дереве нет (документ относится к A7; вместо `STARR`/`STVAR` рекомендуется функция `FIXVAL`).

| файл | сигнатура | назначение | зависимость |
|---|---|---|---|
| `a2eqdsk.f90` | модуль `a2eqdsk`: `EQDSK(coco_number, fileq)` | запись равновесия в G-EQDSK (257×257) в заданном соглашении COCO | — |
| `a2solps.f90` | `A2SOLPS(irad, i_mod)` | запись входного файла для SOLPS (связка ASTRA–SOLPS, один шаг) | SOLPS (внешний) |
| `a2strahl.f90` | модуль `strahl_mod`: `A2STRAHL(tau_start, zneocl, dzneocl, dimpsol)` | связка с кодом переноса примесей STRAHL (R. Dux): запуск бинарника, обмен профилями `Dz/Vz`, источниками, `nimp`, `zavg`, `prad`, `zeff` | STRAHL (закрытый, `$ASTRA_EXT/strahl/sep23/bin`) |
| `a2tglf_elite.f90` | `a2tglf_elite()` | запись геометрии поверхностей для TGLF/ELITE в `tglf/tglf4elite.dat` | — |
| `a2travis.f90` | `A2TRAVIS` | вызов кода ECRH TRAVIS (стелларатор) через его Fortran-API, чтение мощности и тока | TRAVIS (`$(TRAVIS_LIB)`) |
| `ablation.f90` | `ABLATION(trace, pel_prof)` | абляция пеллеты по U-файлам из namelist `&pellet` | — (данные) |
| `ablation_ngs.f90` | `ABLATION_NGS(src, ndep)` | модель абляции пеллеты для стелларатора (Neutral Gas Shielding), трассировка по поверхностям, `&pellet_ngs` | — |
| `alfs.f90` | `alfs(ped_width, a_lfs, dt_tetop, avdte)` | радиус на внешнем обводе (LFS) и нормированный градиент Te в пьедестале | — |
| `dkes_wrapper.f90` | `DKES_INTERFACE_WRAPPER(er_field, no_edge)` | неоклассика стелларатора из таблиц DKES (E. Buglione-Ceresa, 2025) | таблицы DKES, `src/astra/a2dkes.f90` |
| `er_omp.f90` | `er_omp(er_min, er_sep, wexb_lfs, er_lfs, vdia_lfs, bp_lfs)` | профиль Er по величинам на внешней средней плоскости, для TGLF (M. Bergmann, 2026) | использует геометрию SPIDER (по комментарию) |
| `facit.f90` | `facit(Z_imp_in, A_imp_in, N_imp_in, rot_mod, Dz_out, Vz_out)` | неоклассический перенос примесей FACIT (D. Fajardo, 2024) | — |
| `fgauss.f90` | `FGAUSS(YCENTR, YWIDTH, YPROF)` | гауссов профиль осаждения с нормировкой ∫P dV = 1 | — |
| `fgauss_rhopol.f90` | `FGAUSS_RHOPOL(rhop_center, rhop_width, gauss)` | то же по ρ_pol | — |
| `four_arr.f90` | `FOUR_ARR(tim, n_cycle, ARRIN, FREQ, n_harm, arrout)` | разложение Фурье по времени радиального профиля | — |
| `four_scal.f90` | `FOUR_SCAL(tim, n_cycle, VARIN, FREQ, n_harm, VAROUT)` | то же для скаляра | — |
| `foureqc.f90` | `foureqc()` | фурье-моменты магнитных поверхностей в `xpr/fort.four_coef` | — |
| `gnex.f90` | `GNEX()` | поток электронов `GNX` от всех нейтральных источников (Pereverzev 1998) | — |
| `lexthrs.f90` | `lexthrs(lext)` | внешняя индуктивность по Hirshman–Neilson (1986), мкГн | — |
| `mixinq.f90` | `MIXINQ(OPTION, RECOND)` | модель пилообразных колебаний (Кадомцев; Парайль–Переверзев), вариант | — |
| `mixint.f90` | `MIXINT(OPTION, RECOND)` | модель пилообразных колебаний по Кадомцеву (§4 «Sawtooth oscillations») | — |
| `nbi.f90` | `NBI()` | интерфейс к пакету NBI А. Р. Полевого: мощность, импульс, ток от нескольких источников, потери на рипле (§4 «NB heating») | `src/nbi/` (открытый), `exp/nbi/<exp>` |
| `nclass_mod.f90` | модуль `nclass_mod`: `NEOCL4()` | неоклассика NCLASS: χ, D, бутстреп, проводимость (`xd_nc`, `xe_nc`, `jbs_nc`, `cc_nc`, …) | — |
| `neo_ipc.f90` | модуль `a2neo`: `neo_ipc(rho_norm_max)` | вызов NEO через IPC (`xpr/neo`) | NEO, MPI |
| `neut.f90` | `NEUT()` | кинетика нейтралов: `NN`, `TN`, `ALBPL` по `NNCL`, `NNWM`, `ENCL`, `ENWM` (§4 «Gas puff neutrals») | — |
| `pid_control.f90` | `pid_control(Kp_in, Ki_in, Kd_in, error_in, control_signal_in, control_signal_out)` | ПИД-регулятор (используется в `tglf_pid`) | — |
| `qlk_ipc.f90` | модуль `a2qlk`: `qlk_ipc(rho_norm_max)` | вызов QuaLiKiZ через IPC (`xpr/qlki`) | QuaLiKiZ, MPI |
| `qlknn_serial.F90` | `qlknn_serial(chii, chie, e_pflux)` | нейросеть QLKNN в процессе; без `QLKNN_INSTALLED` — заглушка со `stop` | QLKNN |
| `rabbit.F90` | модуль `a2rabbit`: `RABBIT(pNBI_MW, dt_in, pRF_MW, fRF_MHz, nRF_harm, pi_icr, pe_icr)` | NBI-код RABBIT (+ простая модель ИЦР); без `RABBIT_INSTALLED` — заглушка со `stop` | RABBIT (закрытый IPP) |
| `slowdf.f90` | `SLOWDF(pnbi_powers, slow_time)` | торможение быстрых ионов пучка (E. Fable, 2025) | — |
| `slowdn.f90` | `SLOWDN` | время торможения и давление быстрых частиц (альфа), по R. Bilato (E. Fable, 2016) | — |
| `smearr.f90` | `SMEARR(ALFA, f_in, f_out)`, `SMEARRX(…)`, `SMEARR2(…)` | сглаживание профиля минимизацией ∫(α(dU/dx)² + (U−F)²)dx (§4 «Smoothing») | — |
| `stelmetr.f90` | `STELMETR` | отладочная печать баланса стеллараторной метрики (`write(*,*) 'fps', …`) | — |
| `tglf_ipc.f90` | модуль `a2tglf`: `tglf_ipc(rho_norm_max)` | вызов TGLF через IPC (`xpr/tglfi`), до 64 рабочих процессов | TGLF, MPI |
| `tglf_serial.F90` | `tglf_serial(chi_i, chi_e, e_pflux, vimp1, vimp2, i_mflux_as, exchi_as, gamma_as, omega_as)` | TGLF в процессе; без `TGLF_INSTALLED` — заглушка со `stop` | TGLF |
| `torba.F90` | модуль `a2torbeam`: `TORBA(power_MW_in)` | ECRH/ECCD кодом TORBEAM по гиротронам из `&torbeam` → `PEECR`, `CUECR`; без `TORBEAM_INSTALLED` — заглушка со `stop` | TORBEAM (закрытый IPP) |
| `torfpql_mod.F90` | модуль `torfpql_mod`: `toric(pwic, frq, ntor, toll, nmax, debug, prf_abs)` | обёртка TORIC-SSFPQL (ИЦР, R. Bilato); без `TORIC_INSTALLED` — заглушка со `stop` | TORIC |
| `ufr.f90` | `UF1DR(ufnam, tim_in, val)`, `UF2DR(ufnam, tim_in, arr1d)` | чтение 1D/2D U-файла и интерполяция на момент времени | — |

---

## 8. Прочие каталоги

### 8.1. `xpr/` — IPC-программы для TGLF, QuaLiKiZ, NEO

| файл | назначение |
|---|---|
| `tra_interf.c` | общий `main` рабочего процесса: получает из командной строки имя IPC-файла, ключ, номер процесса, размер блока, ID разделяемой памяти и семафора; в зависимости от имени программы вызывает `tglf_interf_`, `qlk_interf_` или `neo_interf_` |
| `tglf_interf.f90` | `tglf_interf(jproc, dims_in, scal_in, prof_in, prof_out)` — запуск TGLF на своём наборе радиусов (`tglf_interface`, `tglf_pkg`) |
| `qlk_interf.f90` | `qlk_interf(…)` — то же для QuaLiKiZ |
| `neo_interf.f90` | `neo_interf(…)` — то же для NEO (`neo_interface`) |

Собираются `exe/Makexpr` в `xpr/tglfi`, `xpr/qlki`, `xpr/neo` (с MPI). Вызывающая сторона — `sbr/*_ipc.f90` и `src/astra/ipc_control.c`. На Ubuntu установщиком отключены.

### 8.2. `preProcess/` — подготовка входных данных

Оба скрипта используют внешний пакет `trview` (`plasma_state`), которого в репозитории нет (IPP).

| файл | назначение |
|---|---|
| `imas2astra.py` | `imas2astra(shot, idsRun, db)`: читает IDS из IMAS (по умолчанию shot 38384, run 6, база `aug`) и пишет вход ASTRA в `$AWD` |
| `sf2astra.py` | `sf2astra(shot, tbeg, tend)`: читает shotfiles AUG (IDA `ne`/`Te`, CEZ `Ti_c`, CEZ/CMZ `vrot`, равновесие EQH), фитирует профили, пишет `exp`/`udb` и namelist RABBIT (`aug_nml` с путями `$ASTRA_EXT/rabbit/...`) |

### 8.3. `pyparse/` (копия в пользовательской области)

`pyparse/iondens_ass.py` — класс `NIAS` с фрагментом Fortran `iondensassign`: `NI(J) = F1(J)+…+F9(J)` («complete AUG»). Подключается основным парсером `pyparse/eqns_init.py` (в корне исходников), если включено `config.checkeqn` — тогда после уравнения `NE` плотность ионов считается суммой `F1…F9`. В корневом `pyparse/` этого файла нет; копия в пользовательской области позволяет пользователю переопределить сборку `NI`.

### 8.4. `exe/` — сборка и запуск

Подробно — [02-build-and-run.md](02-build-and-run.md). Кратко:

| файл | назначение |
|---|---|
| `as_exe` | обёртка bash: определяет `AWD`, платформу, запускает `exe/as_exe.py` |
| `as_exe.py` | разбор ключей (`-m`, `-v`, `-s`, `-e`, `-dev` (машина, по умолчанию `aug`), `-batch`, `-tpause`, `-debug`, `-re`, `-fs`, `-workflow`, `-W`, `-c`, `-nodes`, `-resize`), пересчёт матриц Грина (`greenMatrices.main()`), запись `tmp/astra.nml`, вызов `exe/Build`, затем `json2cdf`. Без аргументов повторяет последний расчёт |
| `Build` | компиляция под модель: `get_platform` + `platform/env.<p>` + `exe/astra_rc[_<compiler>]`, `make -f exe/Makefile`, `exe/Makexpr`, запуск `bin/<m>_gui.exe` или `_batch.exe` |
| `Makefile` | основная сборка: генерация Fortran парсером, библиотеки `libastra`, `libfeqis`, `libnbi`, `libmisc`, `libgraph`, `libsbr`, `libfnc`, линковка с внешними `*_LIB` |
| `Makeastra` | сборка подкаталогов `src/<sub>` в библиотеки (флаг `-DSPIDER` при наличии `SPIDER_LIB`) |
| `Makesub` | сборка `sbr/` и `fnc/` (флаги `-DTORIC_INSTALLED`, `-DTORBEAM_INSTALLED` и т. д. по наличию библиотек) |
| `Makexpr` | сборка IPC-программ `xpr/tglfi`, `qlki`, `neo` |
| `astra_rc` | пути к внешним модулям в `$ASTRA_EXT` (json 9.0.2, TORBEAM nov24, RABBIT may26, TGLF jun25, NEO dec24, QuaLiKiZ/QLKNN nov24, SPIDER latest, TORIC apr26, STELLOPT); `*_LIB` экспортируются только если файл библиотеки существует |
| `astra_rc_gcc`, `astra_rc_ifx` | то же для явного выбора компилятора (`-c gcc` / `-c ifx`) |
| `greenMatrices.py` | расчёт матриц Грина FEQIS: `exp/cnf/<m>_description_in.json` → `_description_full.json` |
| `wr_nml` | вызывает `trview/astra_nml_u.py` через `module load trview` (окружение IPP) |
| `version` | баннер «ASTRA 8.6, January 2026» |
| `README` | что собирается (`mod/`, `lib/`) и что внешние модули выбираются в `exe/astra_rc` |

---

## 9. Наблюдения и открытые вопросы

Дубли и устаревшее:

- `exp/nml/aug34954` и `exp/nml/aug34954_t` побайтно совпадают.
- `exp/nml/steltest2`, `tokstel`, `dtt_exp_R2d19` не имеют парных файлов в `exp/` и ссылаются на отсутствующие данные (`udb/dtt*`, `limiter_walpha.dat`).
- `exp/readme` + `equ/showdata` — материал ASTRA 7, в A8 не запускается (нет `equ/log/showdata`, нет U-файлов).
- `equ/log/stelxes` — в формате «version 7.1.0».
- В `fml/` 12 формул ссылаются на отсутствующие `MIX`/`MIX1`/`SMOFML`, `eli3` — на `elin3R`; `pdtf` помечена устаревшей.
- `equ/tok2vmec` вызывает несуществующую `TOKSTELOUT`.
- `docu/section4.tex` описывает подпрограммы A7 (`MINMAX`, `SETNAV`, `STARR`, `STVAR`, `TSCTRL`, `GNXSRC`), которых в A8 нет.

Открытые вопросы:

- Что делает ядро при `IPEQL=4` без SPIDER: явной остановки в `gs_solver.F90` не найдено, поведение не проверено.
- Запускаются ли `flux_cuas`, `flux_emeq`, `flux_nbi` после комментирования `TORBA`/`RABBIT` и `imep_pw04/08` с `exp/30000_3.4` — не проверено.
- Значение `ZRD3`, `ZRD8` в `exp/aug34954` и суффикса `_08` в `tglf_rot_pred_08` не установлено.
- Мощность NBI в `udb/36982_3400/P36982_3400.NBI_AVG` (10.06 «Вт») выглядит как ошибка единиц.
