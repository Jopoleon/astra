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
| TGLF, NEO (GACODE), QuaLiKiZ, QLKNN | турбулентный и неоклассический транспорт | публичные (GitHub, gitlab.com), нужен MPI | пока не ставятся |

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

`exe/Makexpr` по умолчанию собирает `xpr/tglfi`, `xpr/qlki`, `xpr/neo`. Для них нужны библиотеки из `$ASTRA_EXT` и MPI; `QLK_LIB` в `astra_rc` прописан без проверки существования файла. Без них `Build` завершается с кодом 2 ещё до запуска модели. Upstream-`install.sh` убирает эти цели только в режиме `-safe` и после интерактивных вопросов; установщик убирает их всегда (как ответ `n` в `-safe`).

### 3.7. Тестовая модель требует закрытых кодов

`install.sh` сразу запускает `flux_feqis`/`aug34954`, а эта модель вызывает `TORBA` (ECRH) и `RABBIT` (NBI). Установщик комментирует эти две строки в `equ/flux_feqis` (как `install.sh -safe` при ответе `n`).

### 3.8. Парсер файлов эксперимента

`pyparse/exp_parser.py::parse_u_line` делит строку U-файла по первому **пробелу**. Строки в стиле ASTRA 7 используют табуляцию (`TEX<tab>U-file: 00000.tes<tab>factor: 1.`), а строки-комментарии (`****  ... U-file "00000.tes"`) тоже содержат слово `U-file`. Парсер падал с `IndexError: list index out of range`. Установщик накладывает патч: разбор регулярным выражением, понимающим оба формата; строки без `<VAR> U-file:` пропускаются. Формат A8 (`TE  U-file:34954/TE34954.IDA_AVG: 1.e-3`) разбирается как раньше.

### 3.9. `as_exe` всегда завершается с кодом 0

Даже при ошибке компиляции `exe/as_exe` возвращает 0 (результат `os.system` игнорируется). Об успехе можно судить только по выводу: `>>> ASTRA normal exit >>>` и наличию `ncdf_out/<exp><model>.CDF`. Установщик проверяет именно это.

## 4. Что делает `install_astra8.sh`

| шаг | действие | повторный запуск |
|---|---|---|
| 0 | проверка ОС, подтверждение | — |
| 1 | `apt-get install` пакетов из таблицы раздела 2 | apt пропускает установленное |
| 2 | по HTTPS с LF скачивается **только нужный коммит** (`git fetch --depth 1`, ~90 МБ вместо ~1 ГБ истории; `--full-history` — вся история) | `fetch` + сброс отслеживаемых файлов, затем патчи заново |
| 3 | сборка json-fortran в `$PREFIX/ext/json/9.0.2`; при `--blas mkl` — MKL с PyPI в `$PREFIX/ext/mkl/` | пропускается, если библиотека есть |
| 4 | `env.ubuntu`, `env.mypc`, правка `get_platform`, патч парсера | идемпотентно |
| 5 | копирование `repo_user_area/*` в рабочую копию, отключение xpr-целей и TORBA/RABBIT | `exe/` обновляется, в `equ/`, `exp/`, `udb/` ... дописываются только **отсутствующие** файлы |
| 6 | команда `~/.local/bin/astra8`; если `~/.local/bin` нет в `PATH`, строка добавляется в `~/.bashrc` (стандартный `~/.profile` Ubuntu добавляет его только при следующем входе в систему) | перезаписывается, строка в `.bashrc` не дублируется |
| 7 | `make clean`, полная сборка и batch-расчёт `flux_feqis`/`aug34954`, t = 4–5 с | выполняется заново (`--no-test` — пропустить) |

Параметры:

```
--prefix DIR        корень установки (по умолчанию ~/astra)
--ref REF           ветка, тег или ПОЛНЫЙ хэш коммита ASTRA (по умолчанию проверенный e66f6ba1...;
                    main — свежий код); короткий хэш работает только с --full-history
--full-history      скачать всю историю git ASTRA (~1 ГБ)
--blas openblas|mkl
--no-apt            не ставить пакеты
--no-test           не запускать тестовый расчёт
-y                  без вопросов
```

Журналы: `$PREFIX/install.log` (весь вывод apt/cmake/git), `$PREFIX/test_flux_feqis.log` (сборка и расчёт).

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
- **TGLF / NEO / QuaLiKiZ**: открытые, их можно добавить в установщик (OpenMPI из apt + сборка GACODE/QuaLiKiZ). Пока не сделано.
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
| `FileNotFoundError: .../udb/...` | файл эксперимента ссылается на отсутствующий U-файл | положить U-файл в `~/astra/a8/udb/` |
| окно не появляется, `cannot open display` | нет X-сервера | Windows 11 + WSL2 (WSLg) или рабочий стол Ubuntu; в Windows 10 нужен X-сервер (VcXsrv) и `export DISPLAY=...`; или считать с `-batch` |
| расчёт «молча» не идёт, код возврата 0 | `as_exe` не передаёт ошибки | смотреть вывод: ищите `Error` выше по тексту |
| `cannot fetch '...'` при установке | короткий хэш коммита в `--ref` | дать полный хэш (40 символов) или добавить `--full-history` |
| `astra8: command not found` сразу после установки | терминал открыт до установки | открыть новый терминал |
| Docker: `docker: command not found` в WSL | не включена интеграция Docker Desktop с этой Ubuntu | Docker Desktop → Settings → Resources → WSL integration |
| Docker: `Permission denied` в `~/astra-data` | образ собран с другим UID | пересобрать через `astra-docker.sh` (подставляет ваш UID) или с `--build-arg USER_UID=$(id -u)` |
