# 07. Нагрев, токовая генерация и внешние модули

**О чём глава.** Какие модули нагрева, генерации тока, переноса и примесей ASTRA 8 умеет вызывать, какие из них лежат в дереве `astra-src` целиком, какие представлены только интерфейсом к закрытой или отдельно устанавливаемой библиотеке, как они включаются при сборке (`<MODULE>_LIB` → `-D<MODULE>_INSTALLED`), что происходит, если библиотеки нет (заглушки после коммита `2ec99235`), как устроен механизм параллельного запуска TGLF/NEO/QuaLiKiZ через System V IPC, и что из всего этого реально доступно на «голом» Ubuntu. Отдельно разобраны встроенный NBI, пеллеты и наиболее важные подпрограммы `sbr/`, связанные с нагревом и транспортом. Состояние кода — `astra-src` `67281112` (2026-10-01); при сверке использовалась рабочая установка `~/astra/a8` (только чтение).

Что **не** повторяется здесь:

- синтаксис вызова подпрограмм в модели (`NAME(args)<:dt:tmin:tmax:key;`) и кодогенерация — [03-model-language.md](03-model-language.md);
- сборочный конвейер, `astra_rc`, `Makexpr` в контексте всего запуска — [02-build-and-run.md](02-build-and-run.md);
- SPIDER и FEQIS как решатели равновесия — [06-equilibrium.md](06-equilibrium.md);
- однострочный каталог всех `sbr/fnc/fml` — [09-user-area-catalog.md](09-user-area-catalog.md);
- где оказываются выходные профили нагрева (JSON/CDF) — [08-graphics-and-output.md](08-graphics-and-output.md);
- термины (`PEBM`, `CUECR`, `ZRD`, `CNB` и т. п.) — [glossary.md](glossary.md).

Коды RABBIT, TORBEAM, SPIDER, STRAHL — закрытые коды IPP. Их исходники в дереве отсутствуют, и глава их не описывает и не пытается получить; описано только то, что есть в `astra-src`: интерфейсы, заглушки, данные, скрипты установки.

## Содержание

1. [Сводная таблица модулей](#1-сводная-таблица-модулей)
2. [Общая схема: как модуль попадает в расчёт](#2-общая-схема-как-модуль-попадает-в-расчёт)
3. [Условная компиляция и заглушки (`2ec99235`)](#3-условная-компиляция-и-заглушки-2ec99235)
4. [Общий интерфейс модулей нагрева (по `docu/section4.tex`)](#4-общий-интерфейс-модулей-нагрева-по-docusection4tex)
5. [Встроенный NBI (`src/nbi`, `sbr/nbi.f90`)](#5-встроенный-nbi-srcnbi-sbrnbif90)
6. [RABBIT](#6-rabbit)
7. [TORBEAM (ECRH/ECCD)](#7-torbeam-ecrheccd)
8. [TORIC-SSFPQL (ICRH)](#8-toric-ssfpql-icrh)
9. [LH, внешние источники из данных, гауссовы профили](#9-lh-внешние-источники-из-данных-гауссовы-профили)
10. [SPIDER](#10-spider)
11. [STRAHL (примеси)](#11-strahl-примеси)
12. [Транспортные модели: TGLF, NEO, QuaLiKiZ, QLKNN](#12-транспортные-модели-tglf-neo-qualikiz-qlknn)
13. [Механизм IPC: `ipc_control.c`, `xpr/`, `Makexpr`](#13-механизм-ipc-ipc_controlc-xpr-makexpr)
14. [Пеллеты: `ABLATION`, `ABLATION_NGS`, модель `flux_pel`](#14-пеллеты-ablation-ablation_ngs-модель-flux_pel)
15. [Прочие важные подпрограммы `sbr/`](#15-прочие-важные-подпрограммы-sbr)
16. [Установка модулей: `astra_rc`, `modules_install.sh`, `modules_inst_loc.sh`](#16-установка-модулей-astra_rc-modules_installsh-modules_inst_locsh)
17. [Документация A7 против реальности A8](#17-документация-a7-против-реальности-a8)
18. [Что делать на «голом» Ubuntu](#18-что-делать-на-голом-ubuntu)
19. [Где смотреть в исходниках](#19-где-смотреть-в-исходниках)

---

## 1. Сводная таблица модулей

«Открыт» — исходник физики лежит в дереве и компилируется всегда. «Интерфейс» — в дереве только Fortran-обёртка, сама физика во внешней библиотеке или программе. «Ubuntu» — что нужно сделать на обычной машине без доступа к GitLab IPP.

| Модуль | Назначение | Лицензия / статус | Что в дереве | Как включается | Без библиотеки | Ubuntu |
|---|---|---|---|---|---|---|
| NBI (Polevoi) | нагрев/ток/частицы от пучков | открыт | `src/nbi/*` + `repo_user_area/sbr/nbi.f90` | всегда (`lib/libnbi.a`) | — | работает из коробки; нужен файл `exp/nbi/<exp>` или `exp/nbi/<machine>` |
| RABBIT | быстрый Монте-Карло-подобный NBI (+ простая ICRF) | закрыт (IPP, `gitlab.mpcdf.mpg.de:markusw/rabbit`) | `repo_user_area/sbr/rabbit.F90`, данные `rabbit/limiters`, `rabbit/tables_*` | `RABBIT_LIB` → `-DRABBIT_INSTALLED` | заглушка печатает сообщение и делает `stop` | недоступен; вызов `RABBIT` в модели надо закомментировать |
| TORBEAM | ECRH/ECCD, трассировка пучков | закрыт (IPP) | `repo_user_area/sbr/torba.F90` | `TORBEAM_LIB` → `-DTORBEAM_INSTALLED` | `stop` | недоступен; `TORBA` закомментировать; ECRH задавать вручную (§9) |
| TORIC-SSFPQL | ICRH, полноволновой + квазилинейный ФП | внешний (R. Bilato, IPP), исходников нет | `repo_user_area/sbr/torfpql_mod.F90`, `python/eqdsk_to_TORIC_gs.py` | `TORIC_LIB` → `-DTORIC_INSTALLED` | `stop` | недоступен |
| SPIDER | свободно-граничное равновесие | закрыт (IPP) | вызовы в `src/astra/gs_solver.F90`, `stepup.F90` | `SPIDER_LIB` → `-DSPIDER` (`exe/Makeastra`) | вызов просто не компилируется (§10) | недоступен; использовать FEQIS (см. [06](06-equilibrium.md)) |
| STRAHL | транспорт и излучение примесей | закрыт (IPP, `rld/strahl`) | `repo_user_area/sbr/a2strahl.f90` (всегда компилируется), атомные данные `strahl/` | внешний исполняемый файл `$ASTRA_EXT/strahl/sep23/bin/strahl` | `system()` вернёт ошибку, файлы результатов не появятся | недоступен |
| TGLF | квазилинейный турбулентный транспорт | открыт (GACODE, MIT-подобная лицензия GA) | `sbr/tglf_ipc.f90`, `sbr/tglf_serial.F90`, `xpr/tglf_interf.f90` | `TGLF_LIB` + `Makexpr` (IPC-программа `xpr/tglfi`) | IPC: программы нет → зависание/ошибка; serial: `stop` | можно собрать самостоятельно (§12), но в коде есть дефекты |
| NEO | неоклассика (drift-kinetic) | открыт (GACODE) | `sbr/neo_ipc.f90`, `xpr/neo_interf.f90` | `NEO_LIB` + `Makexpr` (`xpr/neo`, нужен MPI) | IPC-программы нет | можно собрать самостоятельно |
| QuaLiKiZ | квазилинейный гирокинетический транспорт | открыт (gitlab.com/qualikiz-group) | `sbr/qlk_ipc.f90`, `xpr/qlk_interf.f90` | `QLK_LIB` + `Makexpr` (`xpr/qlki`, нужен MPI) | IPC-программы нет | можно собрать самостоятельно |
| QLKNN | нейросетевой суррогат QuaLiKiZ | открыт (gitlab.com) | `sbr/qlknn_serial.F90` | `QLKNN_LIB` → `-DQLKNN_INSTALLED` | `stop` | можно собрать самостоятельно |
| NCLASS / NEOCL4 | неоклассика | открыт, исходник в дереве | `sbr/nclass_mod.f90` | всегда | — | работает |
| FACIT | неоклассический транспорт примесей | открыт, в дереве | `sbr/facit.f90` | всегда | — | работает (не проверено прогоном) |
| TRAVIS | ECRH для стеллараторов | внешний | `sbr/a2travis.f90` | `TRAVIS_LIB` нигде не определён | — | недоступен |
| DKES | неоклассика стеллараторов | внешний | `sbr/dkes_wrapper.f90`, `postProcess/plotdkesrun.py` | внешняя программа | — | вне задач токамака |
| SOLPS | граничная плазма | внешний | `sbr/a2solps.f90` (обмен файлами) | — | — | вне задач |
| Пеллеты | абляция/депозиция | открыт | `sbr/ablation.f90`, `sbr/ablation_ngs.f90` | всегда | — | работает при наличии файлов профиля абляции |
| json-fortran | ввод/вывод JSON | открыт | — (внешняя библиотека) | `JSON_HOME` в `astra_rc` | сборка невозможна | обязательна; ставится установщиком |

Ключевой практический вывод: **на чистом Ubuntu из модулей нагрева работает только встроенный NBI**. ECRH/ICRH надо либо задавать аналитически (гауссовы профили, §9), либо брать из внешних расчётов как экспериментальные данные (`CARnX`).

---

## 2. Общая схема: как модуль попадает в расчёт

1. **Вызов в модели.** Модуль вызывается в разделе уравнений модели как подпрограмма, например `NBI<:0.05:::;`, `TORBA<:0.05:::;`, `tglf_ipc(CF16)<:0.02:::;`. Разбор полей `:dt:tmin:tmax:key` и префиксов `<<`, `<`, `>` выполняет `parse_sbr` в `pyparse/parse_as.py` (≈441–500):
   - `<<` → `locsbr=-2`;
   - `<` → `locsbr=-1`: вызывается в `DETVAR` до равновесия;
   - `>` (а также всегда для `MIXINT`, `MIXEXT`, `TSCTRL`) → `locsbr=1`: после шага (postep);
   - без префикса → `locsbr=0`: внутри итераций шага.
   Подробно — [03-model-language.md](03-model-language.md).
2. **Генерируемый код** всегда подключает модули-обёртки, независимо от того, вызывает ли их модель: `pyparse/const_text.py` (≈200–205, 895–904) вставляет `use a2rabbit, only: rabbit`, `use a2torbeam, only: torba`, `use torfpql_mod, only: toric`, `use a2tglf`, `use a2qlk`, `use a2neo`, `use strahl_mod`, `use nclass_mod`, а в `DETVAR` безусловно вызываются `tglf_alloc`, `qlk_alloc`, `neo_alloc`. Поэтому обёртки должны компилироваться всегда — отсюда и заглушки (§3).
3. **Входные данные** модуль берёт напрямую из модулей Fortran ядра (`scalars`, `status`, `parameters_a2equil` и т. д.), а свои настройки — из:
   - namelist-файла `exp/nml/<exp>` (если нет — `exp/nml/<machine>`; `src/astra/read_input.f90` ≈145–149);
   - для NBI — файла `exp/nbi/<exp>` или `exp/nbi/<machine>` (§5);
   - U-файлов в `udb/<shot>/` (мощности, углы гиротронов и т. п.).
4. **Выход** записывается в стандартные массивы ядра (`PEBM`, `PIBM`, `CUBM`, `PEECR`, `CUECR`, `PEICR`, `PIICR`, `PEFW`, `PIFW`, `CUFW`, ...), которые модель затем явно подставляет в источники: `PE=POH+PEBM+PEECR-...`, `CD=CUBM+CUECR` (см. `repo_user_area/equ/flux_feqis`).

Полный список «слотов» нагрева в `astra_variables.json`:

| Канал | Мощность в электроны | В ионы | Ток | Скаляры |
|---|---|---|---|---|
| NBI | `PEBM` | `PIBM` | `CUBM` | `QNBI` |
| ECR | `PEECR` | — | `CUECR` | `QECR`, `FECR` |
| ICR | `PEICR` | `PIICR` | `CUICR` | `QICR`, `FICR` |
| FW (быстрая волна) | `PEFW` | `PIFW` | `CUFW` | `QFW`, `FFW` |
| LH | `PELH` | — | `CULH` | `QLH`, `FLH` |

Плюс NBI-специфичные `NIBM`, `NNBM1..3`, `SNEBM`, `SNIBM1..3`, `SNNBM`, `SCUBM`, `PBLON`, `PBPER`, `PFAST`, `PELON`.

---

## 3. Условная компиляция и заглушки (`2ec99235`)

### 3.1. Флаги

Флаги задаются в `repo_user_area/exe/Makesub` по наличию переменной окружения с путём к библиотеке:

| Переменная | Флаг | Файл | Строки `Makesub` |
|---|---|---|---|
| `TORIC_LIB` | `-DTORIC_INSTALLED` | `sbr/torfpql_mod.F90` | 15–19, 92 |
| `TORBEAM_LIB` | `-DTORBEAM_INSTALLED` | `sbr/torba.F90` | 22–26, 86 |
| `RABBIT_LIB` | `-DRABBIT_INSTALLED` | `sbr/rabbit.F90` | 29–33, 89 |
| `QLKNN_LIB` | `-DQLKNN_INSTALLED` | `sbr/qlknn_serial.F90` | 36–40, 83 |
| `TGLF_LIB` | `$(TGLF_FLAG)` — **не определён** | `sbr/tglf_serial.F90` | 80 |
| `SPIDER_LIB` | `-DSPIDER` | `src/astra/gs_solver.F90`, `stepup.F90` | `exe/Makeastra` 13–17, 52 |

Сами `*_LIB` экспортирует `repo_user_area/exe/astra_rc` (и варианты `astra_rc_*`): для каждого модуля задаётся путь в `$ASTRA_EXT/<module>/<version>/...`, и если файл библиотеки не существует, переменная сбрасывается (`unset`). Таким образом флаг появляется только при реально существующем файле.

**Дефект:** правило для `tglf_serial.o` использует `$(TGLF_FLAG)`, но переменная `TGLF_FLAG` не определяется ни в `Makesub`, ни в `Makefile`, ни в `astra_rc` (`grep -rn TGLF_FLAG repo_user_area/exe` находит только строку 80). Следовательно, `tglf_serial.F90` **всегда компилируется как заглушка**, даже если TGLF установлен (вывод по коду; прогоном с установленным TGLF не проверено). IPC-вариант `tglf_ipc` от этого не зависит.

### 3.2. Коммит `2ec99235`

Коммит Giovanni Tardini от 2026-05-12, «Failsafe for modules SPIDER, TGLF, QLKNN, RABBIT, TORBEAM, TORIC in case of missing installation»:

- переименовал `rabbit`, `torba`, `tglf_serial`, `qlknn_serial`, `torfpql_mod` в `.F90` (чтобы работал препроцессор);
- обернул реальный код в `#ifdef <MODULE>_INSTALLED ... #else <заглушка> #endif`;
- поправил `Makesub`, `astra_rc*`, `Build`, `Makexpr`, `ipc_control.c`.

Все заглушки устроены одинаково: печатают строку вида `RABBIT not installed, check RABBIT_LIB path in exe/astra_rc` и вызывают `stop`:

| Подпрограмма | Заглушка | Строки |
|---|---|---|
| `RABBIT` | `repo_user_area/sbr/rabbit.F90` | 433–442 |
| `TORBA` | `repo_user_area/sbr/torba.F90` | 437–445 |
| `TORIC` | `repo_user_area/sbr/torfpql_mod.F90` | 603–612 |
| `QLKNN_SERIAL` | `repo_user_area/sbr/qlknn_serial.F90` | ветка `#else` |
| `TGLF_SERIAL` | `repo_user_area/sbr/tglf_serial.F90` | ветка `#else` |

Следствие: «failsafe» означает «сборка проходит», а не «расчёт продолжается». Модель, вызывающая `TORBA` или `RABBIT`, на машине без этих библиотек **соберётся, но остановится при первом вызове**. Поэтому `install.sh -safe` и наш установщик комментируют эти вызовы в `equ/flux_feqis` (`sed -i "s#RABBIT#\!RABBIT#g"` и аналогично для `TORBA`), см. [02-build-and-run.md](02-build-and-run.md) и `installer/docs/DETAILED.md`.

**Несоответствие интерфейса TORIC.** Реальная `toric(pwic, frq, ntor, toll, nmax, debug, prf_abs)` имеет 7 аргументов, заглушка — 6 (без `prf_abs`). Для модели `flux_toric` (`TORIC(pwic=1.0d0,nmax=2,debug=0)`) это не важно, но вызов с `prf_abs` без TORIC не скомпилируется.

---

## 4. Общий интерфейс модулей нагрева (по `docu/section4.tex`)

Раздел «Interfaces to additional heating and CD modules» (`docu/section4.tex` ≈2842–3208) — документация ещё ASTRA 6/7, но общий принцип в A8 сохранён.

**Входные данные**, которые модуль получает из ASTRA: геометрия (`RTOR`, `AB`/`ABC`, `AMETR`, `SHIF`, `ELON`, `TRIA`), состав плазмы (`AMJ`, `ZMJ`, `ZEF`, `ZIM1..3`, `NE`, `NI`, `NHYDR`, `NDEUT`, `NIZ1..3`), температуры `TE`, `TI`, ток и поле (`IPL`, `BTOR`, `MU`, `CU`).

**Выходные данные**: плотности мощности в электроны/ионы `P_e^H`, `P_i^H` и плотность тока увлечения `j_CD`.

Документ подчёркивает, что все дополнительные пакеты, кроме NBI, предоставляются «by special request», а NBI включён всегда. В A8 это правило буквально отражено в сборке: NBI — `lib/libnbi.a`, остальные — внешние библиотеки по `*_LIB`.

---

## 5. Встроенный NBI (`src/nbi`, `sbr/nbi.f90`)

Автор физики — A. R. Polevoi. Это единственный модуль нагрева с полным открытым исходником в дереве.

### 5.1. Файлы

| Файл | Роль |
|---|---|
| `repo_user_area/sbr/nbi.f90` (80 строк) | точка входа `NBI` для модели: `set_input` → `NBINJ(trim(NBFILE), ...)` → `get_output` |
| `src/nbi/nb_injection.f90` | `NBINJ`: чтение конфигурации, расчёт депозиции, замедления, источников |
| `src/nbi/nbibce.f90` и др. | сечения, вспомогательные процедуры; `nbibce.f90` использует `dat/fij.dat`, `dat/srsfi.dat` |
| `src/nbi/ReadMe.txt` | заметки об изменениях A8 (переназначение выходов, `CNB4`, `CNBI4`) |

### 5.2. Конфигурационный файл

Путь вычисляется при чтении входных данных: `src/astra/read_input.f90:82` вызывает `inquire_fname('nbi', exp_file, machine, NBFILE)` (`src/astra/parse_utils.f90` ≈177–195):

1. `exp/nbi/<exp_file>`, если существует;
2. иначе `exp/nbi/<machine в нижнем регистре>`;
3. иначе `NBFILE='***'`.

Если `NBFILE` начинается с `*`, `nbi.f90` печатает `>>> NBI Error >>> Configuration file not found` и останавливает расчёт. Интерактивного диалога (как в A7) нет.

Формат строк конфигурации — поля по 12 символов (`docu/section4.tex`, табл. 4.31; в A8 `NBINJ` читает 20 полей через `STREAD`):

| № | Поле | Единицы / значение |
|---|---|---|
| 1 | мощность источника | МВт (обычно переопределяется `ZRDn`, см. ниже) |
| 2 | доля counter-инжекции | 0…1 |
| 3 | масса частиц пучка | а.е.м. |
| 4 | заряд | — |
| 5 | максимальная энергия | кэВ |
| 6–8 | доли мощности на E, E/2, E/3 | — |
| 9 | усреднение по орбитам | 0/1/2 |
| 10 | число «тонких» пучков | — |
| 11 | вертикальный сдвиг | м |
| 12–13 | `Rmax`, `Rmin` | м |
| 14 | `tan α` | — |
| 15 | аспект пятна | — |
| 16–19 | `CVER1/2`, `CHOR1/2` — параметры распределений `NBFHZ`/`NBFRY` из `sbr/nbuser.f` | в A8 не используются (`nbuser.f` в дереве нет) |
| 20(21) | резерв | — |

Пример — `repo_user_area/exp/nbi/aug34954` (8 источников AUG). Мощности по источникам в A8 задаются через `ZRD1..ZRD8`; в `exp/aug34954` заданы только `ZRD3 2.5` и `ZRD8 2.5`, т. е. работают два источника по 2.5 МВт.

### 5.3. Управляющие параметры

Из табл. 4.30 `docu/section4.tex` и кода `nb_injection.f90`:

| Параметр | Смысл в документации | Фактически в A8 |
|---|---|---|
| `CNB1` | число источников | `n_nbi<0` → `n_nbi+1`; `n_nbi==0` → `stop` |
| `CNB2` | потери на перезарядку (1 — явно) | используется |
| `CNB3` | прореживание сетки `N1=(NA1-1)/CNB3+1` | **игнорируется**: `dn_rho` внутри сначала полагается 1, затем при `fp_flag>0` вычисляется из ларморовского радиуса |
| `CNB4` | 1 — стационарный, 2 — нестационарный ФП | по `ReadMe.txt`: `CNB4=0` — источник для TORIC в `dat/toric.nbi`; `CNB4=2` — функция распределения в `dat/nbion.dat`, если `dat/nbi.dat` содержит `1` |
| `CNBI1..4` | доп. параметры | `CNBI4` управляет термоядерными реакциями пучок–плазма (`ReadMe.txt`) |

Модели перезарядки 1/2/3 (документация): `PBICX=-0.0024*SNNBM*TI`.

### 5.4. Выходы

`CUBM`, `CUFI`, `NIBM`, `NNBM1..3`, `PBEAM`, `PBLON`, `PBPER`, `PEBM`, `PIBM`, `SCUBM`, `SNEBM`, `SNIBM1..3`, `SNNBM` (табл. 4.32). В A8 часть выходов переиспользована (`nb_injection.f90` ≈232–249, `ReadMe.txt`): `NNBM2`, `NNBM3`, `SNIBM1..3` содержат источники нейтронов/трития (`stnbdp`, `sdnbdp1/2`, `sdnbtp`), а не плотности горячих нейтралов, как в A7. Изменения «August 2026» добавили `PBCX`, `SBICX`, `SRSTH`.

### 5.5. Модель `flux_nbi`

`repo_user_area/equ/flux_nbi` (строки ≈44–60):

```
CNB1=8; CNB2=1; CNB3=2; CNB4=1;
NBI<:0.05:::;
TORBA<:0.05:::;
EQDSK(7)<:0.05:::;
```

То есть модель **одновременно** вызывает встроенный NBI и TORBEAM. На машине без TORBEAM вызов `TORBA` надо закомментировать, иначе расчёт остановится в заглушке. Регрессия `regressions/aug34954flux_nbi.CDF`: максимум `PIBM` ≈0.76 МВт/м³ (против 1.34 в `flux_feqis` с RABBIT), `PEECR` присутствует (TORBEAM). Прогон `flux_nbi` с закомментированным `TORBA` на gfortran в рамках этой главы **не проверялся**.

---

## 6. RABBIT

Модуль `a2rabbit` в `repo_user_area/sbr/rabbit.F90`.

**Реальная ветка** (`#ifdef RABBIT_INSTALLED`): `RABBIT(pNBI_MW, dt_in, pRF_MW, fRF_MHz, nRF_harm, pi_icr, pe_icr)` использует `mod_rabbit_lib` (`rabbit_lib_init`, шаг по времени, `get_icrh_depo` и др.).

- Настройки — namelist'ы `&rabbit_beam_geo`, `&partmix`, `&nbi_par`, `&physics`, `&species` из `exp/nml/<exp>` (пример — `repo_user_area/exp/nml/aug34954`).
- Контур лимитера — `rabbit/limiters/limiter_aug.dat` (в дереве есть также `COMU`, `d3d`, `dtt`, `future`, `jet`, `tcv`, `vns`).
- Таблицы — `rabbit/tables_highRes/` (есть также `tables_ITERDEMO/`).
- Требует `NHYDR`/`NDEUT`/`NTRIT`; сечения для примеси есть только для `Zimp` = 4, 5, 6, 7, 10, 28, иначе принудительно 6 (углерод).
- Выходы: `PEBM`, `PIBM`, `NIBM`, `CUBM`, `SNEBM`, `SCUBM`, `PBLON`, `PBPER`, `NRATE`; плюс простая модель ICRF через `pRF` (`pi_icr`, `pe_icr`).

**В дереве** — только данные: `rabbit/` (≈13 МБ: `limiters/`, `tables_*`), исходников нет.

**Заглушка** (строки 433–442) — сообщение и `stop`.

---

## 7. TORBEAM (ECRH/ECCD)

Модуль `a2torbeam` в `repo_user_area/sbr/torba.F90`.

- Вызов: `TORBA(power_MW_in)` (в моделях обычно без аргументов: `TORBA<:0.05:::;`).
- Namelist `&torbeam` в `exp/nml/<exp>`: `n_gyro`, `pecr_file`, `theta_file`, `phi_file` (U-файлы мощности и углов, напр. `udb/34954/P34954.ECH_AVG`, `THE34954.ECH_AVG`, `PHI34954.ECH_AVG`), `beam_on`, `nmod`, `freq_n`, `theta_n`, `phi_n`, `RR_n`, `ZZ_n`, `xryyb`, `xrzzb`, `xwyyb`, `xwzzb`.
- До 30 гиротронов (`n_gy_max=30`); A7-документация (`docu/additional_packages_add_EF.tex`) описывает `TORBA(VAR1,VAR2)`, 8 гиротронов и файл в `exp/ecr/` — **устарело**.
- Равновесие передаётся на прямоугольной сетке 64×64 (`ctr2rz_b`).
- Выходы: `PEECR`, `CUECR`.
- Заглушка (437–445) — `stop`.
- Библиотека по `astra_rc`: `libtorbeamB.so` + `-lirc` (Intel runtime), версия `nov24`.

---

## 8. TORIC-SSFPQL (ICRH)

Обёртка R. Bilato, `repo_user_area/sbr/torfpql_mod.F90`.

Последовательность реальной ветки `toric(pwic, frq, ntor, toll, nmax, debug, prf_abs)`:

1. создаёт каталог `AWD/toric/torfpql_<дата>_<время>`;
2. пишет `torica.inp` из шаблона `exp/icrh/torica.inp` с подстановкой namelist `&toric` (`freqcy`, `hf_pwgoal`, `accur_goal`, `max_iter`, `nphi`) из `exp/nml/<exp>`;
3. пишет EQDSK через `a2eqdsk::eqdsk` (COCOS 11);
4. вызывает `python3 ../../python/eqdsk_to_TORIC_gs.py toric.eqdsk --rho_max 0.99` (скрипт 3018 строк; нужны `numpy`, `scipy`, `matplotlib`, `scikit-image`);
5. запускает TORIC, читает поглощение;
6. выходы: `PIICR`, `PEICR`, `PIFW`, `PEFW`, `CUFW`, `CUICR` (плотность передаётся в единицах 10¹³ см⁻³).

Пути TORIC в `astra_rc` жёстко указывают на MPCDF (`FFTW_LIBDIR=/mpcdf/soft/...`, `TORIC .../apr26/intel/2025.1`) — вне IPP не переносится без правки. Модель `repo_user_area/equ/flux_toric` вызывает `TORBA`, `RABBIT(pRF...)` и `TORIC(pwic=1.0d0,nmax=2,debug=0)`.

---

## 9. LH, внешние источники из данных, гауссовы профили

**LH-модуля в A8 нет.** Есть только слоты `PELH`, `CULH`, `QLH`, `FLH`. Для LH, а на Ubuntu и для ECRH/ICRH, остаются три пути:

1. **Аналитические профили.** `sbr/fgauss.f90`: подпрограмма `FGAUSS(YCENTR, YWIDTH, YPROF)` заполняет `YPROF(j)=exp(-((ρ_j-YCENTR·ROC)/(YWIDTH·ROC))²)` и нормирует так, что `∫YPROF dV = 1`. Пример из комментария в файле: `FGAUSS(CF1, CF2, CAR1):; PE=...+CF3*CAR1;` — тогда `CF3` есть полная мощность в МВт. Есть вариант по ρ_pol (`sbr/fgauss_rhopol.f90`, в `equ/tglf` закомментирован вызов `FGAUSS_RHOPOL(0.5d0, 0.3d0, CAR54)>:.2:2.:5.:;`) и сглаживание `sbr/smearr.f90`. Для произвольного профиля мощность нормируется интегралом `VINT`, как в `flux_pel`: `CAR53 = CF16*CAR51/VINT(CAR51B)`.
2. **Источники из внешнего расчёта как экспериментальные данные.** Профиль, посчитанный другим кодом (TRANSP/NUBEAM, внешний трассировщик ЭЦ), кладётся в U-файл и подключается в `exp/<exp>` как `CARn U-file:...`, затем в модели подставляется `CARnX`. Так устроена модель `equ/imep_pw04` (строки ≈40–48): `PEBM=car20x; PIBM=car21x; SNEBM=CAR18X; NIBM=car19x; PEECR=0.0;`.
3. **Встроенный NBI** для пучков (§5).

---

## 10. SPIDER

Свободно-граничный решатель равновесия (закрыт). Подробнее как решатель — [06-equilibrium.md](06-equilibrium.md); здесь — только аспект «внешнего модуля».

- Флаг `-DSPIDER` ставится в `repo_user_area/exe/Makeastra` (13–17, 52) при наличии `SPIDER_LIB` (`libspider.a`, версия `latest`).
- `src/astra/gs_solver.F90` ≈700–836 (`A_EQUIL`): параметры SPIDER — под `#ifdef SPIDER`; если `equil_solver==101`, вызывается FEQIS (`feqis_main`), иначе `spider_run` — **только** под `#ifdef SPIDER`.
- `src/astra/stepup.F90:229`: `SPIDUPDATE` под тем же `#ifdef`.
- Namelist `&spider` — в `exp/nml/<exp>`.

Заглушки со `stop` для SPIDER нет. Модели с `IPEQL=4` (`flux_spider`, `flux_neo_tglf`, `astrahl_simple_msp`) без SPIDER, по-видимому, будут работать **без пересчёта равновесия** (вызов просто отсутствует). Это вывод по коду, прогоном не проверено. Такой режим опаснее явной остановки: расчёт пойдёт, но равновесие будет начальным.

---

## 11. STRAHL (примеси)

Интерфейс открыт и компилируется всегда (`repo_user_area/sbr/a2strahl.f90`, модуль `strahl_mod`), сам STRAHL — внешняя закрытая программа.

- Вызов: `A2STRAHL(tau_start, zneocl, dzneocl, dimpsol)`, например в `equ/astrahl_simple_msp`: `A2STRAHL(0.d0,1.d6,1.d6,CSOL1)>::::;`.
- Рабочий каталог: копирует `AWD/strahl/` (`atomdat`, `*.atomdat`, `pec_files`, `pstrahl`, `strahl.control`) в `AWD/strahl/<exp>_<equ>/`.
- Запускает через `system()` внешний `$ASTRA_EXT/strahl/sep23/bin/strahl a q` (или `$ASTRA_EXT/strahl/bin/strahl`, если `ASTRA_EXT==AWD/astra_modules`), затем `result_to_astra`; разбирает `results.txt`.
- Namelist `&strahl_par` в `exp/nml/<exp>`.
- Массивы результата: `zeff_strahl`, `prad_tot_strahl`, `nmain_strahl`, `prad_main_strahl`, `prad_strahl`, `nimp_strahl`, `zavg_strahl`, `nesrc_strahl`, `Dneo_strahl`, `Vneo_strahl`, `Dz_in_strahl`, `Vz_in_strahl`, `rrates_in_strahl`. Группа `strahl` описана в `astra_variables.json` (13 переменных), но в JSON-выход **не пишется** (см. [08-graphics-and-output.md](08-graphics-and-output.md)).
- В дереве `strahl/` (≈8 МБ) — только атомные данные (`Ar`, `B_`, `Be`, `C_`, `Cu`, `Fe`, `H_`, `He`, `Kr`, `Li`, `N_`, `Ne`, `Ni`, `O_`, `W_`, `Xe` `.atomdat`, `atomdat/newdat`, `pec_files`, `pstrahl`, `strahl.control`).

Явной проверки кода возврата `system()` и наличия бинарника перед запуском в обёртке не видно; поведение при отсутствии STRAHL — ошибка чтения `results.txt` (по коду, не проверено).

---

## 12. Транспортные модели: TGLF, NEO, QuaLiKiZ, QLKNN

Это не нагрев, но они входят в ту же категорию внешних модулей и используют тот же механизм сборки. Для TGLF, NEO и QuaLiKiZ основной путь — **параллельные IPC-программы** (§13), для TGLF и QLKNN есть также последовательные варианты.

### 12.1. Обёртки

| Обёртка | `ipcId` | Входов на точку | Выходов | Программа | Замечания |
|---|---|---|---|---|---|
| `sbr/tglf_ipc.f90` (`a2tglf`) | 0 | 47 (`n_inputs`) | 15 (`n_arr_out`) | `xpr/tglfi` | `n_dims=7`, `n_scalars=8`; геометрия `geom_flag=1` (ELITE, опция 3) |
| `sbr/neo_ipc.f90` (`a2neo`) | 1 | 33 | — | `xpr/neo` | MPI |
| `sbr/qlk_ipc.f90` (`a2qlk`) | 2 | 39 | — | `xpr/qlki` | MPI |
| `sbr/tglf_serial.F90` | — | 9 аргументов | — | линкуется напрямую | всегда заглушка (дефект `TGLF_FLAG`, §3.1) |
| `sbr/qlknn_serial.F90` | — | `qlknn_serial(CAR21,CAR22,CAR24)` | — | линкуется напрямую | `QLKNN_INSTALLED` |

Общие для IPC-обёрток параметры: `nrho_m=64` радиальных точек, `nworkers=64` процессов. Профили интерполируются на 64-точечную сетку до `rho_norm_max*ROC`.

**Жёсткое требование:** обёртка останавливается (`stop`), если `nworkers > MAX_NWORKERS` или `mod(nrho_m, nworkers) ≠ 0`. `MAX_NWORKERS` берётся из окружения, иначе `omp_get_num_procs()/2` (`sbr/tglf_ipc.f90` ≈88–106; сообщение `>>> ERROR nworkers= 64 exceeds #physical CPUs=...`). `repo_user_area/exe/Build` при отправке через `sbatch` пишет в скрипт `export MAX_NWORKERS=128` (строка 99), а при любом локальном запуске (GUI или batch) **безусловно** экспортирует `MAX_NWORKERS = nproc / threads_per_core` (строки 109–118), перекрывая значение, заданное пользователем. На рабочей машине разработчика (24 логических ядра, 2 потока на ядро) получается 12 < 64, и **все три IPC-обёртки остановятся**. Чтобы запустить их на обычном ПК, нужно править `Build` (или `nworkers` в обёртках, сохраняя `mod(64, nworkers)=0`); 64 процесса на 12 ядрах будут медленными, но допустимы.

### 12.2. Выходы

- TGLF: `tglf_out%chi_i`, `chi_e`, `e_pflux`, `ion_pflux`, `mom_flux`, `equipart`, `gamma`, `omega` (строки `prof_out` 1..7+). В модели `equ/tglf`: `tglf_ipc(CF16)<:0.02:::;`, затем `CAR21="tglf_out%chi_i(j)"` и т. д.
- NEO: `chii`, `chie`, `jbs`, `elec_pflux`, `vippd`, `vittd`, `vippi1`, `vitti1`.
- Программы-обработчики: `xpr/tglf_interf.f90` → `tglf_run` (рабочий путь `tglf/`), `xpr/qlk_interf.f90` → `qualikiz`, `xpr/neo_interf.f90` → `neo_init_serial('neo/')` + `neo_run`.

### 12.3. Модели, использующие транспортные модули

| Модель | Вызовы |
|---|---|
| `tglf` | `tglf_ipc(CF16)<:0.02:::;` |
| `flux_tglf` | `TGLF_IPC<:0.1:::;` |
| `tglf_pid` | `PID_CONTROL` + `tglf_ipc` |
| `tglf_serial`, `flux_tglf_serial` | `tglf_serial(...)` (9 аргументов) |
| `tglf_rot_pred_08` | TGLF + вращение |
| `qlk` | `qlk_ipc(CF16)` |
| `qlknn` | `qlknn_serial(CAR21,CAR22,CAR24)` |
| `flux_neo` | `NEO_IPC` |
| `flux_neo_tglf` | `IPEQL=4` (SPIDER) + `TGLF_IPC` + `NEO_IPC` |

Большинство этих моделей дополнительно вызывают `TORBA<:.05:::;` и `RABBIT<:.05:::;`.

### 12.4. Библиотеки по `astra_rc`

TGLF `jun25` (`libtglf.a`); NEO `dec24` (`neo_lib.a nclass_lib.a UMFPACK_lib.a math_lib.a expro_lib.a geo_lib.a`; для NEO переменная **не сбрасывается** при отсутствии файлов); QuaLiKiZ `nov24` (`libQLK-intel-release-default.a` / `libQLK-gcc-release-default.a`); QLKNN `nov24`. Источники — GACODE (`github.com/gafusion/gacode`, ветки TGLF и NEO), QuaLiKiZ и QLKNN-fortran с gitlab.com (см. §16).

---

## 13. Механизм IPC: `ipc_control.c`, `xpr/`, `Makexpr`

### 13.1. Сборка программ

`repo_user_area/exe/Makexpr` собирает из `xpr/tra_interf.c` (с `-Dtglf`, `-Dneo`, `-Dqlk`) и соответствующего `*_interf.f90` три исполняемых файла: `xpr/tglfi`, `xpr/neo`, `xpr/qlki`. Для `neo` и `qlki` нужен MPI-компилятор (`MPIFC`). Создаются рабочие каталоги `AWD/tglf`, `AWD/qualikiz`, `AWD/neo`. Статические библиотеки MKL прописаны жёстко. `exe/Build` вызывает `Makexpr` после основной сборки. Наш установщик удаляет эти цели, т. к. библиотек нет.

### 13.2. Сторона ASTRA: `src/astra/ipc_control.c`

| Функция (Fortran-имя) | Что делает |
|---|---|
| `initialise_ipc_` | `ftok` по пути исполняемого файла ASTRA и PID; `key2 = key + ipcId*3`; создаёт `Nsub+1` семафоров и 3 сегмента разделяемой памяти (`dims`, `vars`, `arrs`); пишет `tmp/<exp><equ>-<ipcId>.ipc`; для каждого рабочего процесса через `system()` запускает в фоне `AWD/xpr/<name> <ipc_file> <key> <j> <N_CHUNK> <N_ARR_OUT> <shm ids> &`; ждёт семафор 0 |
| `fill_int_shm_`, `fill_dbl_shm_` | копирует входные данные в общую память |
| `unlock_sbp_` | разрешает рабочим считать |
| `wait4all_` | ожидание через `semtimedop` с таймаутом 0.1 с (на macOS — эмуляция) |
| `sbp2astra_` | читает id выходного сегмента рабочего из ipc-файла и забирает результат |

### 13.3. Сторона рабочего процесса: `xpr/tra_interf.c`

`main` подключается к сегментам ASTRA, создаёт собственный сегмент (`N_CHUNK*N_ARR_OUT` double), дописывает в ipc-файл строку `J_PROC PID ShmId`, затем в цикле: `semop` на семафоре 0 (+1, «готов»), ожидание своего семафора, вызов `tglf_interf_`/`neo_interf_`/`qlk_interf_(J_PROC, dims, scal, prof_in, prof_out)`, сброс семафора. Процесс завершается, когда семафор удалён.

### 13.4. Очистка

`exe/as_exe.py` после расчёта читает `tmp/<expequ>-<0..2>.ipc`, убивает перечисленные PID и удаляет сегменты памяти `ipcrm -m`. Семафоры при этом **не удаляются** (по коду). При аварийном завершении могут остаться и процессы, и сегменты; в GUI есть клавиша `J` (`ipcs -s`, `ipcs -m`) для диагностики ([08-graphics-and-output.md](08-graphics-and-output.md)).

### 13.5. Старая документация

`docu/IPC.tex` (эпоха ASTRA 6.1) описывает синтаксис `NEUT & : 0.1 : 2.2 : 3.3;`, шаблон `IPC/template.c`, цепочку `toexample/to_example/ot_example/otexample`, файл `AWD/tmp/astra.ipc`, `Nsbp+1` семафоров и `Nsbp+2` сегментов. `docu/icp_add_EF.tex` описывает `xpr/subp(j1,j2)&::::;` с разбиением по радиусу и приоритеты `<`/`>`. В A8 этого синтаксиса нет: `pyparse` не знает `&`, а модель вызывает `tglf_ipc` как обычную подпрограмму, которая сама запускает рабочих. Механизм семафоров и общей памяти концептуально тот же.

---

## 14. Пеллеты: `ABLATION`, `ABLATION_NGS`, модель `flux_pel`

### 14.1. `sbr/ablation.f90`

`ABLATION(trace, pel_prof)`:

- namelist `&pellet` (`rho_abl_file`, `time_abl_file`, `mass`) в `exp/nml/<exp>`;
- чтение профиля абляции из U-файлов (`ufheader`, `ufrd`, `uf1dr`);
- возвращает профиль источника частиц `pel_prof`.

В `repo_user_area/exp/nml/aug34954` namelist `&pellet` **отсутствует** (там есть `rabbit_beam_geo`, `partmix`, `nbi_par`, `species`, `output_settings`, `physics`, `numerics`, `torbeam`, `spider`, `equil_settings`, `strahl_par`, `toric`). Поэтому `flux_pel` с `aug34954` из коробки, вероятно, упадёт на чтении namelist (не проверено). Для своего разряда нужен свой `&pellet` и U-файлы абляции.

### 14.2. Модель `flux_pel`

`ABLATION(CF16, CAR51)>::::;` — вызов после шага; нормировка `CAR53 = CF16*CAR51/VINT(CAR51B);`, затем `CAR53` подставляется в источник частиц.

### 14.3. `sbr/ablation_ngs.f90`

`ABLATION_NGS(src, ndep)` — модель NGS/MCL23 для стеллараторов, читает `PEL_CHORD_FILE`; сопутствующий скрипт `vmec/python/pellet.py`. Для токамака Глобус-М2 не предназначена.

---

## 15. Прочие важные подпрограммы `sbr/`

Только связанные с нагревом, источниками и транспортом; полный каталог — [09-user-area-catalog.md](09-user-area-catalog.md).

| Файл | Что делает | Внешние зависимости |
|---|---|---|
| `nclass_mod.f90` | `NEOCL4` + полный исходник NCLASS (бутстреп, неоклассические коэффициенты); используется в `imep_pw04/08` (`jbs_nc`) | нет |
| `facit.f90` | FACIT, неоклассический транспорт примесей (D. Fajardo, 2024) | нет |
| `mixint.f90`, `mixinq.f90` | модель пилообразных колебаний Кадомцева (Переверзев); вызывается после шага | нет |
| `neut.f90` | нейтралы (`NEUT<:0.05:::;` в `flux_nbi`) | нет |
| `slowdf.f90`, `slowdn.f90` | замедление быстрых ионов (Fable) | нет |
| `alfs.f90` | альфа-частицы | нет |
| `lexthrs.f90` | порог L–H | нет |
| `pid_control.f90` | ПИД-регулятор (`tglf_pid`) | нет |
| `er_omp.f90` | радиальное поле Er для TGLF (M. Bergmann, 2026) | нет |
| `a2tglf_elite.f90` | геометрия ELITE для TGLF | TGLF |
| `a2eqdsk.f90` | запись EQDSK; `EQDSK(coco)` пишет `ncdf_out/<exp><equ><время>.eqdsk` | нет |
| `fgauss.f90`, `fgauss_rhopol.f90`, `smearr.f90` | гауссовы и сглаженные профили депозиции | нет |
| `ufr.f90` | чтение U-файлов `UF1DR`/`UF2DR` | нет |
| `gnex.f90` | вспомогательный ввод | нет |
| `four_*.f90`, `foureqc.f90`, `stelmetr.f90` | фурье-описание поверхностей, стеллараторная метрика | нет |
| `dkes_wrapper.f90` | DKES (стелларатор) | внешний DKES |
| `a2travis.f90` | TRAVIS ECRH (стелларатор), вызывает внешний `set_nTZ_profs_f77` | `TRAVIS_LIB` нигде не задаётся |
| `a2solps.f90` | обмен файлами с SOLPS | SOLPS |

---

## 16. Установка модулей: `astra_rc`, `modules_install.sh`, `modules_inst_loc.sh`

### 16.1. `repo_user_area/exe/astra_rc`

Версии на 2026-10-01: `JSON_HOME=$ASTRA_EXT/json/9.0.2`, TORBEAM `nov24`, RABBIT `may26`, TGLF `jun25`, NEO `dec24`, QuaLiKiZ `nov24`, QLKNN `nov24`, SPIDER `latest`, TORIC `apr26/intel/2025.1`. Для `darwin` используются каталоги без версии. Общая раскладка — `$ASTRA_EXT/<module>/<version>/lib/...`. Подробно — [02-build-and-run.md](02-build-and-run.md).

### 16.2. `modules_install.sh` (кластер IPP)

Ставит Intel oneAPI, cmake, json-fortran 9.0.2, цепочку NetCDF, затем модули из git:

| Модуль | Источник |
|---|---|
| RABBIT | `git@gitlab.mpcdf.mpg.de:markusw/rabbit.git` (закрытый) |
| TORBEAM | `git@gitlab.mpcdf.mpg.de:ipp-aug/torbeam` (закрытый) |
| SPIDER | `git@gitlab.mpcdf.mpg.de:git/spider` (закрытый) |
| STRAHL | `git@gitlab.mpcdf.mpg.de:rld/strahl.git` (закрытый) |
| QuaLiKiZ | `https://gitlab.com/qualikiz-group/QuaLiKiZ.git` (открытый, MPI) |
| QLKNN-fortran | gitlab.com + 4 репозитория namelist'ов (открытый) |
| GACODE (TGLF, NEO) | `git@github.com:gafusion/gacode.git` (открытый) |

### 16.3. `modules_inst_loc.sh` (локальная машина)

Вариант для linux (apt) / darwin (brew). TORBEAM — через `ssh://gerrit.ipp.mpg.de:29418/ipp/e1/ecrh/libtorbeam`, STRAHL — `azito/strahl`, GACODE — через conda-окружение `CONDA_OMPI_GNU`. Известные несоответствия (по коду, не проверено прогоном):

- модули ставятся в каталоги без версии `$ASTRA_EXT/<module>`, тогда как `astra_rc` для не-darwin Linux ищет версионированные пути — библиотеки не будут найдены без правки `astra_rc`;
- скрипт подключает `platform/env.linux`, которого нет (есть `aug columbia cz darwin docker freia gway hgw iter lac mit ncku omega perlmutter puhti rat2 sevilla tohtori tok tok_gcc tok_ifx w7x`; `env.darwin` задаёт `ASTRA_EXT=$AWD/astra_modules`).

### 16.4. `requirements.txt`

Компилятор Intel > 2018 или gfortran; MKL, MPI, OpenMP; NetCDF «если RABBIT»; Python ≥ 3, `scipy ≥ 0.19.1`. Наш установщик ставит NetCDF «на будущее» и закрытые модули не ставит (`installer/docs/DETAILED.md`).

---

## 17. Документация A7 против реальности A8

| Тема | `docu/*.tex` (A6/A7) | A8 (код) |
|---|---|---|
| NBI-конфигурация | `exp/<exp>.nbi` или `exp/aug.nbi`, интерактивный диалог | `exp/nbi/<exp>` или `exp/nbi/<machine>`, диалога нет, без файла — `stop` |
| `CNB3` | прореживание радиальной сетки | игнорируется |
| `NNBM2/3`, `SNIBM1..3` | горячие нейтралы, источники ионов по энергиям | частично переназначены под нейтронные/тритиевые источники |
| TORBEAM | `TORBA(VAR1,VAR2)`, 8 гиротронов, `exp/ecr/` | `TORBA(power)`, до 30 гиротронов, namelist `&torbeam` в `exp/nml/` |
| TORIC | `a2toric.f` | `torfpql_mod.F90`, TORIC-SSFPQL, Python-конвертер равновесия |
| STRAHL | `A2STRAHL.f` | `a2strahl.f90`, внешний бинарник `$ASTRA_EXT/strahl/sep23/bin/strahl` |
| IPC | `NAME & :...;`, `IPC/template.c`, `AWD/tmp/astra.ipc` | обычный вызов `tglf_ipc`, `tmp/<exp><equ>-<id>.ipc`, программы `xpr/tglfi|neo|qlki` |
| Отсутствие пакета | пакет просто не поставлялся | обёртка компилируется, вызов → `stop` (кроме SPIDER) |

`docu/additional_packages_add_EF.tex` также описывает рестарт `RESTARTA/IRESTA`; в A8 рестарт выполняется из JSON (`as_exe -re N`), см. [08-graphics-and-output.md](08-graphics-and-output.md).

---

## 18. Что делать на «голом» Ubuntu

1. **Ничего не подключать «на всякий случай».** `astra_rc` сам сбросит отсутствующие `*_LIB`; сборка пройдёт с заглушками.
2. **Закомментировать в модели** все вызовы `TORBA`, `RABBIT`, `TORIC`, `tglf_serial`, `qlknn_serial` (`!` в начале строки). Это делает наш установщик для `flux_feqis`; для других моделей — вручную.
3. **Не использовать `IPEQL=4`** (SPIDER) — модель пойдёт без пересчёта равновесия (предположительно); переключиться на FEQIS ([06-equilibrium.md](06-equilibrium.md)).
4. **NBI** — встроенный: свой файл `exp/nbi/<exp>` (или `exp/nbi/<machine>`) по образцу `exp/nbi/aug34954`, мощности через `ZRDn`, `CNB1` — число источников.
5. **ECRH/ICRH/LH** — гауссовы профили (`FGAUSS`) или профили из внешнего расчёта через `CARnX` (§9).
6. **TGLF/NEO/QuaLiKiZ** — возможны после самостоятельной сборки GACODE/QuaLiKiZ и правки `astra_rc` (пути) и `Makexpr` (MKL → OpenBLAS/LAPACK), плюс `MAX_NWORKERS≥64` через правку `exe/Build` (§12.1). Для `tglf_serial` дополнительно нужно определить `TGLF_FLAG=-DTGLF_INSTALLED` в `Makesub`. Ничего из этого здесь не проверялось.
7. **Регрессии с RABBIT/TORBEAM** (`flux_feqis`, `flux_nbi`, `tglf_pid` и др.) на такой машине численно не совпадут с эталонами: в эталонах есть `PEECR`, `PIBM`, `PIICR` от закрытых кодов ([08-graphics-and-output.md](08-graphics-and-output.md), раздел о регрессиях).

---

## 19. Где смотреть в исходниках

| Что | Где |
|---|---|
| Флаги модулей | `repo_user_area/exe/Makesub` (15–40, 79–101), `repo_user_area/exe/Makeastra` (13–17, 52) |
| Пути и версии библиотек | `repo_user_area/exe/astra_rc`, `astra_rc_*` |
| IPC-программы | `repo_user_area/exe/Makexpr`, `repo_user_area/xpr/{tra_interf.c, tglf_interf.f90, qlk_interf.f90, neo_interf.f90}` |
| IPC на стороне ASTRA | `src/astra/ipc_control.c` |
| Очистка IPC | `repo_user_area/exe/as_exe.py` (после `os.system(... Build ...)`) |
| `MAX_NWORKERS` | `repo_user_area/exe/Build`, `repo_user_area/sbr/{tglf_ipc,neo_ipc,qlk_ipc}.f90` |
| Встроенный NBI | `src/nbi/nb_injection.f90` (≈232–249 — переназначение выходов), `src/nbi/ReadMe.txt`, `repo_user_area/sbr/nbi.f90` |
| Путь к NBI-файлу и namelist'ам | `src/astra/read_input.f90:82`, ≈145–149; `src/astra/parse_utils.f90` ≈177–195 (`inquire_fname`) |
| Пример конфигурации | `repo_user_area/exp/nbi/aug34954`, `repo_user_area/exp/nml/aug34954`, `repo_user_area/exp/aug34954` |
| RABBIT | `repo_user_area/sbr/rabbit.F90` (заглушка 433–442), `rabbit/` |
| TORBEAM | `repo_user_area/sbr/torba.F90` (заглушка 437–445) |
| TORIC | `repo_user_area/sbr/torfpql_mod.F90` (заглушка 603–612), `python/eqdsk_to_TORIC_gs.py` |
| SPIDER | `src/astra/gs_solver.F90` ≈700–836, `src/astra/stepup.F90:229` |
| STRAHL | `repo_user_area/sbr/a2strahl.f90`, `strahl/` |
| Пеллеты | `repo_user_area/sbr/ablation.f90`, `ablation_ngs.f90`, `repo_user_area/equ/flux_pel` |
| Обязательные `use` в сгенерированном коде | `pyparse/const_text.py` ≈200–205, 895–904 |
| Разбор вызова подпрограммы | `pyparse/parse_as.py` ≈441–500 (`parse_sbr`) |
| Список слотов нагрева | `astra_variables.json` (группы `profiles`, `variables`) |
| Документация A6/A7 | `docu/section4.tex` ≈2842–3208, `docu/additional_packages_add_EF.tex`, `docu/IPC.tex`, `docu/icp_add_EF.tex` |
| Скрипты установки модулей | `modules_install.sh`, `modules_inst_loc.sh`, `requirements.txt` |
| Коммит с заглушками | `git show 2ec99235` |
