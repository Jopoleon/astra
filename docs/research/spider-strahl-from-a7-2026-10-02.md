# SPIDER и STRAHL из ASTRA 7 в ASTRA 8 (Ubuntu, gfortran)

Дата: 2026-10-02. ASTRA 8: `e66f6ba1` (рабочая копия `~/astra/a8`, скопирована в отдельный каталог; оригинал не трогался). gfortran 13.3, OpenBLAS, Ubuntu 24.04 (WSL2).

Исходники: каталог лаборатории `nikita_input/` (в `.gitignore`, **закрытые коды IPP, в репозиторий и в публикации не попадают**): полная ASTRA 7.02 (`a7/ASTRA_7.02`, в т. ч. `EQUIL/SPIDER/SRC`, `MAstraN/for`) и пакет `ASTRA_STRAHL/STRAHL`. Никита (2026-10): «Spider не должен был меняться с 7й астры, можно выдрать оттуда (папка EQUIL). Старый strahl тоже можно попытаться взять из 7й астры».

Здесь нет фрагментов закрытого кода; описаны интерфейсы, правки (sed-однострочники) и результаты.

## Итог

| код | результат |
|---|---|
| SPIDER (равновесие, `IPEQL=4`, фиксированная граница) | **работает.** Собирается gfortran в `libspider.a` + наш переходник (≈120 строк). Модели `flux_spider`, `astrahl_simple_msp`, синтетика «Глобуса» доходят до `ASTRA normal exit`; совпадение с FEQIS 0.1–4 %, с эталоном IPP (настоящий SPIDER `jun26`, Intel) ≤ 1.5 % по метрике |
| SPIDER со свободной границей (`ITFBE`/`IFBEY>0`) | **не проверено.** Нужны файлы машины в формате SPIDER A7 (`exp/equ/coil.dat`, `durs.dat`, …); у «Глобуса» в A7 их нет. Сигнатура `SPIDUPDATE` в A8 и A7 разная (см. ниже) |
| STRAHL (примеси, `A2STRAHL`) | **работает.** A7-STRAHL собирается как отдельная программа; наш конвертер параметров (2 отличия формата) и наш `result_to_astra` (Python, netCDF → текст). `astrahl_simple_msp` / `aug34954`, t=4–5 с: совпадение с эталоном IPP ≤ 0.5 % (`NIZ1`, `NIZ2`, `ZEF`, `PRAD`) |
| STRAHL: NEOART, источник электронов | **не поддержано:** `Dneo`, `Vneo`, `nesrc` возвращаются нулями (в выходе этой версии STRAHL их нет); NEOART в A7-версии привязан к общим блокам ASTRA 7 |

В установщик добавлены ключи `--a7-src`, `--spider-src`, `--strahl-src` (секция `5b`, разметка `# >>> a7modules`). Переходники — `installer/a7modules/` (наш код, без кода IPP).

## 1. SPIDER

### 1.1. Что ждёт ASTRA 8

- `exe/astra_rc`: `SPIDER_HOME=$ASTRA_EXT/spider/jun26`, `SPIDER_LIB=$SPIDER_HOME/lib/libspider.a`, `SPIDER_INC=-I$SPIDER_HOME/inc`; если файла нет — `SPIDER_LIB` сбрасывается.
- `exe/Makeastra`: при заданном `SPIDER_LIB` все `*.F90` компилируются с `-DSPIDER $(SPIDER_INC)`; `exe/Makefile` линкует `$(SPIDER_LIB)` после `-lmisc`.
- `src/astra/gs_solver.F90` (`A_EQUIL`, стр. 714–836): `use spider_params, only: type_parameters`, поля `dt time neql ntheta kpr k_grid epsro enels key_plc key_dmf key_out k_fixfree no_circuit_eq key_start nstep`; вызов `spider_run(ncoils, ucoils, equil_in, equil_out, parameters_equil)` с `type_equilibrium` из `src/misc/imas_ids.f90`.
- `src/astra/stepup.F90:229`: `SPIDUPDATE(machine, CCOIL, time, n_coils)` — только при `IFBEY>=1` и `ICIRCQ>0`.

### 1.2. Что есть в A7

`EQUIL/SPIDER/SRC` (Makefile под ifort, цель `esp.a`): `spider_run(ncoils, vcoils, equil_in, equil_out, params)` в `interfac.f`, модули `schemas` (тип `type_equilibrium`, «минимальные CPO» ITM) и `parameters` (`type_parameters`). Вызывалось из `EQUIL/COUPLING_SCHEME/a_spider.f`.

Внимание: копия `nikita_input/ASTRA_7.02/EQUIL/SPIDER/SRC` неполная (нет `interfac.f` и др.); собирать из полной ASTRA 7 (`a7/ASTRA_7.02`).

### 1.3. Различия и как они закрыты

| различие | решение |
|---|---|
| A8 `type_equilibrium` (`imas_ids`) — надмножество A7 `schemas` (добавлены `rho_tor_norm`, `squareness`, `darea`, `jphi`, `fdia2d`, `psiaxis`, `psiext`, …); `teta2d` → `theta2d` | свой `schemas.f90`: `use imas_ids` + константы `RKIND`, `TWOPI`; `sed 's/%teta2d/%theta2d/g' interfac.f` |
| A8 ждёт модуль `spider_params` с полями `ntheta`, `no_circuit_eq`; в A7 — модуль `parameters` с `nteta` | свой `spider_params.f90`: `type, extends(<A7 type_parameters>) :: type_parameters` + `ntheta`, `no_circuit_eq` |
| `equil_out` в `A_EQUIL` — `intent(out)`: указатели обнуляются, A7 заполняет только свои поля; дальше `json_rw.f90` пишет `rho_tor_norm`, `theta2d`, `darea`, `jphi`, `fdia2d` → иначе обращение к пустым указателям | A7-вход переименован (`sed` в `interfac.f`: `subroutine spider_run(` → `spider_run_a7(`); наш `spider_run` в `spider_a8_adapter.f90` передаёт `nteta=ntheta`, вызывает A7 и дозаполняет недостающие поля (нули/нормированная сетка) |
| `PSIEXT`: A7 брал его вызовом `psib_ext`/`f_psib_ext` после расчёта, A8 читает `equil_out%global_param%psiext` | переходник вызывает `psib_ext` (`k_grid=0`) / `f_psib_ext` (`k_grid=1`) при `k_fixfree=1` |
| `SPIDUPDATE`: в A8 4 аргумента, в A7 без аргументов (данные из файлов) | линкуется (внешняя процедура без интерфейса, лишние аргументы игнорируются), но токи катушек из `CCOIL` в A7-SPIDER не попадают. Касается только свободной границы с уравнениями цепей — не проверено |
| `no_circuit_eq` (статическая свободная граница A8) | в A7 аналога нет, игнорируется |
| `k_grid=1` (адаптивная сетка; так в `namelist_spider.txt` лаборатории) | A8 принудительно ставит `k_grids = 0` после чтения `&spider` (`gs_solver.F90`); адаптивная сетка из A8 недоступна без правки A8 |

Компиляция: все объекты `esp.a` из A7 Makefile + переходник, `gfortran -O2 -fPIC -w -std=legacy -fallow-argument-mismatch`. Других правок исходников A7 не нужно; ошибок компиляции нет. Неразрешённых внешних символов и пересечений с библиотеками A8 нет (проверено `nm`). `imas_ids.f90` из A8 компилируется рядом только ради `.mod` (в `libspider.a` не входит).

SPIDER в фиксированной границе пишет служебные файлы в `exp/equ/` (`outp.wr`, `q.wr`, `tab_bnd.wr`, …); без этого каталога расчёт падает (`Cannot open file 'exp/equ/outp.wr'`) — установщик создаёт его (`mkdir -p exp/equ`).

### 1.4. Проверка

Сборка и запуск в копии установки (`ASTRA_EXT` → копия `ext`), `libastra` пересобран с `-DSPIDER`:

```bash
# модель = flux_spider без TORBA/RABBIT (закрытые), т. е. flux_feqis с IPEQL = 4
sed -e 's/^TORBA/!TORBA/' -e 's/^RABBIT/!RABBIT/' equ/flux_spider > equ/spider_test
cp equ/log/flux_spider equ/log/spider_test
astra8 -m spider_test -v aug34954 -s 4 -e 5 -batch
astra8 -m spider_test -v gm2synth -s 0.1 -e 0.5 -dev gm2 -batch
astra8 -m astrahl_simple_msp -v aug34954 -s 4 -e 5 -batch     # SPIDER + STRAHL
```

Все — `>>> ASTRA normal exit >>>` + CDF; с `EQDSK(7)` пишутся `.eqdsk`.

Сравнение (максимум |Δ| / максимум |эталона|, последний момент):

| случай | эталон | `MU` | `CU` | `ELON` | `TRIA` | `VOLUM` | `G11` | `G33` | `SHIF` |
|---|---|---|---|---|---|---|---|---|---|
| `aug34954`, t=5 | FEQIS (`flux_feqis`), наш | 3.6 % | 3.7 % | 0.17 % | 3.3 % | 0.02 % | 2.3 % | 0.08 % | 0.7 % |
| `gm2synth` (A=1.5, κ=1.9), t=0.5 | FEQIS, наш | 0.7 % | 0.5 % | 0.25 % | 1.6 % | 0.06 % | 0.5 % | 0.07 % | 0.3 % |
| `astrahl_simple_msp`/`aug34954`, t=5 | IPP `regressions/aug34954astrahl_simple_msp.CDF` | 0 (задан) | 1.5 % | 0.5 % | 1.4 % | 0.09 % | — | 0.25 % | 4.7 % (на оси 0.052 / 0.054 м) |

SPIDER быстрее FEQIS: `aug34954` 4–5 с — 6 с CPU против 30 с; `gm2synth` — 5 с против 13 с.

## 2. STRAHL

### 2.1. Что ждёт ASTRA 8

`sbr/a2strahl.f90`: STRAHL — **внешняя программа**. В каталоге `strahl/<exp>_<model>/` (копия `$AWD/strahl`: атомные данные A8) пишутся `param_files/sparams.dat`, `nete/grid_11111.1`, `nete/pp11111.1`, затем `cd <dir> && $ASTRA_EXT/strahl/sep23/bin/strahl a q` и `.../result_to_astra <первый элемент> > results.txt`. `results.txt` читается списком: `Nr`, затем блоки «метка — массив»: `rho_pol`, `Zeff`, `Prad_tot`, `n_main`, `Prad_<sp>`, `Prad_main`, `n_<sp>`, `Zavg_<sp>`, `Dneo_<sp>`, `Vneo_<sp>`, `nesrc_<sp>` (единицы СИ: Вт/м³, м⁻³).

### 2.2. Что есть в A7

`ASTRA_STRAHL/STRAHL/str` — STRAHL (Р. Дукс) встроен в ASTRA 7 как подпрограмма: `start_strahl` (аргументы командной строки закомментированы, `astracall=.true.`), связь с ASTRA через общие блоки (`sbr/as7_strahl_4imp.f`, `astraout`, `initstrahl`, …). Часть файлов включает `for/parameter.inc`, `for/const.inc`, `for/status.inc` из ядра ASTRA 7 (`MAstraN/for`). При `astracall=.false.` работает классический путь: чтение файлов, выход — netCDF `result/<El>strahl_result.dat` на каждый элемент. Программы `result_to_astra` в пакете нет.

### 2.3. Что сделано

| шаг | как |
|---|---|
| отдельная программа | наш `strahl_main.f`: разбор аргументов (`a`, `q`, …), вызов `strahl(...)` с `astracall=.false.`; заглушки `astraout`, `astraset`, `astra_tend`, `initstrahl`, `initgrid2` (путь ASTRA 7, не вызывается) |
| компиляция | объекты из `str/makefile` + драйвер; `-I. -I.. -I<A7>/MAstraN` (только ради `for/*.inc` в `neoarchive.f`); `netcdf.inc` из системы; `-lnetcdff -lnetcdf` (`libnetcdff-dev` уже ставится установщиком) |
| правка исходника | одна: `sed 's/neoclass\.eq\.\.false\./neoclass.eqv..false./' read_grid.f` (сравнение логических через `.eq.` — расширение ifort) |
| формат `sparams.dat` | A8 пишет формат STRAHL `sep23`. Отличия для A7-читателя: (1) строка `cv finite vol=1, finite diff=0` + одно число → у A7 `max_step max_diff` (берём `1000 -1.` — как во всех примерах пакета); (2) целые поля, записанные A8 как вещественные (флаг «rate from file», «divertor puff», «prompt redep») → gfortran не читает `0.000` в integer. Конвертер `sparams_a8_to_a7.py` (идемпотентный) запускается обёрткой `strahl` перед `strahl_a7` |
| атомные данные | берутся из `$AWD/strahl` (поставка A8); A7-читатель их понимает (лишний столбец `nbicx` игнорируется) |
| `result_to_astra` | наш `result_to_astra.py`: последний момент из netCDF всех элементов (`impurity_density` по зарядам, `impurity_radiation`: индекс `nion+2` — сумма, `nion-1` — тормозное основной плазмы), `n_main` и `Zeff` — из квазинейтральности по всем примесям, `Zavg = Σz·n_z / Σn_z`. `Dneo`, `Vneo`, `nesrc` = 0 |

### 2.4. Проверка

`astrahl_simple_msp` / `aug34954` (в `exp/nml/aug34954` `nimp_touse = 2`: B и W; Ar модели остаётся нулём), t=4–5 с, 101 шаг, 101 вызов STRAHL, 25 с CPU, `ASTRA normal exit`. Сравнение с эталоном IPP (STRAHL `sep23` + SPIDER `jun26`, Intel), t=5:

| `NIZ1` (B) | `NIZ2` (W) | `ZIM1` | `ZIM2` | `ZEF` | `PRAD` | `CAR51` (Prad main) | `CAR53` | `CAR54` | `F1`/`F2` (встроенные уравнения) |
|---|---|---|---|---|---|---|---|---|---|
| 0.49 % | 0.31 % | 0.01 % | 0.09 % | 0.15 % | 0.33 % | 0.04 % | 0.5 % | 0.34 % | 0.47 % / 0.52 % |

## 3. Установщик

```bash
bash installer/install_astra8.sh -y \
    --a7-src  /path/to/ASTRA_7.02 \
    --strahl-src /path/to/ASTRA_STRAHL/STRAHL
```

- `--a7-src DIR` — SPIDER из `DIR/EQUIL/SPIDER/SRC`; `DIR/MAstraN/for` нужен STRAHL. `--spider-src DIR` — явно каталог SPIDER. `--strahl-src DIR` — требует `--a7-src`.
- Исходники копируются в `$PREFIX/build/{spider,strahl}` и собираются там; результат: `$PREFIX/ext/spider/<версия из exe/astra_rc>/{lib/libspider.a,inc/*.mod}`, `$PREFIX/ext/strahl/<версия из sbr/a2strahl.f90>/bin/{strahl,strahl_a7,result_to_astra,*.py}`. Версии читаются из файлов ASTRA, а не зашиты.
- Полная пересборка ASTRA при каждой установке уже есть (`make clean`), поэтому `-DSPIDER` применяется сам.
- В `equ/flux_spider` комментируются `TORBA`/`RABBIT`. Тесты (если не `--no-test`): `flux_spider`/`aug34954` 4–5 с; при SPIDER+STRAHL — `astrahl_simple_msp`/`aug34954` 4–4.2 с. Логи `$PREFIX/test_<model>.log`.
- Без этих ключей поведение установщика не меняется.

## 4. Ограничения и ловушки

- **Свободная граница SPIDER не проверена**; для «Глобуса» нет файлов машины SPIDER, а `SPIDUPDATE` не получает токи из A8.
- **Адаптивная сетка SPIDER (`k_grid=1`) в A8 выключена** самим A8.
- **NEOART в STRAHL** не работает в этой связке (нули). Источник электронов от ионизации примесей — нули.
- `nimp_touse` в `exp/nml/<exp>` должен совпадать с числом примесей модели; лишние остаются нулями.
- **Ловушка A8, не связанная со SPIDER:** имя файла `.eqdsk` собирается в `character(len=120)` из `AWD` (`sbr/a2eqdsk.f90:90`); при длинном пути установки (> ~70 символов) `EQDSK` падает с `End of record`. `AWD` зашивается в `src/tmp/ininam.f90` при сборке: после переноса каталога нужно пересобрать модель (`touch equ/<модель>`).
- `as_exe` всегда возвращает 0 — успех определяется по `ASTRA normal exit` и CDF.

## 5. Что лежит где

| файл | что | в git |
|---|---|---|
| `installer/a7modules/spider/schemas.f90` | замена A7-модуля `schemas` через `imas_ids` | да (наш код) |
| `installer/a7modules/spider/spider_params.f90` | тип параметров A8 как расширение A7-типа | да |
| `installer/a7modules/spider/spider_a8_adapter.f90` | `spider_run` для A8 | да |
| `installer/a7modules/strahl/strahl_main.f` | драйвер-программа + заглушки | да |
| `installer/a7modules/strahl/strahl.sh`, `sparams_a8_to_a7.py` | обёртка `strahl`, конвертер параметров | да |
| `installer/a7modules/strahl/result_to_astra.{sh,py}` | netCDF → `results.txt` | да |
| исходники SPIDER/STRAHL | только у пользователя (`--a7-src`, `--strahl-src`) | **нет** |

## 6. Проверка установщика в чистых контейнерах

Исходники A7 подключались только томом на чтение (в образ и в репозиторий не попадают):

```bash
docker run --rm -v "$PWD/installer:/inst:ro" -v "$PWD/nikita_input:/a7:ro" ubuntu:24.04 \
  bash -c 'bash /inst/install_astra8.sh -y --a7-src /a7/a7/ASTRA_7.02 --strahl-src /a7/ASTRA_STRAHL/STRAHL; echo EXIT=$?'
```

| образ | gfortran | SPIDER | STRAHL | `flux_feqis` | `flux_spider` 4–5 с | `astrahl_simple_msp` 4–4.2 с | EXIT |
|---|---|---|---|---|---|---|---|
| `ubuntu:24.04` | 13.3 | собран | собран | normal exit | normal exit | normal exit, 21 вызов STRAHL | 0 |
| `ubuntu:22.04` | 11 | собран | собран | normal exit | normal exit | normal exit | 0 |
