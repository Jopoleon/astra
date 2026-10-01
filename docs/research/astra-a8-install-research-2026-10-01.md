# ASTRA 8: почему не ставится на «обычный» Ubuntu — исследование

Дата: 2026-10-01. Upstream проверен на коммите `e66f6ba1` (`main`, последняя активность 2026-10-01).

## Запрос

Лаборатория (Никита Жильцов, 2026-09-29/30): ASTRA 8 из `https://gitlab.mpcdf.mpg.de/git/astra` не удаётся нормально установить и запустить. Нужна **инструкция по развёртыванию на чистой ОС** (любой Ubuntu).

Критерии «работает» (со слов лаборатории):

- аналог hello world: модель `showdata` на данных `readme`; в A7 это `astra readme showdata`, в A8 — `exe/as_exe -m showdata -v readme`;
- пример из репозитория: `exe/as_exe -m flux_feqis -v aug34954 -s 4 -e 5`;
- ожидаемое поведение и картинки есть для старой ASTRA 7.02, которая у них работает (скриншоты от 2026-09-30, X11-окно ASTRA с профилями).

Файл `install A8.txt` с тем, что Никита уже делал, в рабочую область не попал — **его нужно запросить**.

## Ссылки

| Ссылка | Статус |
|---|---|
| `https://gitlab.mpcdf.mpg.de/git/astra` | публичный, клонируется по HTTPS без аккаунта; ~1.2 ГБ, ~1800 коммитов, LGPL-2.1; мейнтейнеры G. Tardini, E. Fable (IPP) |
| `README.md` | советует `git clone git@...` (SSH): без аккаунта MPCDF это не работает, нужен HTTPS: `git clone https://gitlab.mpcdf.mpg.de/git/astra.git a8` |
| RABBIT `gitlab.mpcdf.mpg.de:markusw/rabbit`, SPIDER `git/spider`, TORBEAM `ipp-aug/torbeam` / gerrit IPP, STRAHL `azito/strahl` | **закрытые** (API отдаёт 404 анонимно), только SSH для аккаунтов IPP |
| QuaLiKiz, QLKNN (`gitlab.com/qualikiz-group`), json-fortran, NetCDF, CMake | публичные |
| GACODE (TGLF/NEO) `git@github.com:gafusion/gacode.git` | публичный, но в скрипте задан SSH-URL: без SSH-ключа GitHub клон упадёт, нужен HTTPS |

## Как устроена установка

```
install.sh
 ├─ get_platform          → имя платформы по hostname/пути (aug, tok, iter, ..., darwin, иначе "mypc")
 ├─ source platform/env.<platform>  → FC, CC, MPIFC, FC_FLAGS, MKL_LIBDIR, ASTRA_EXT
 ├─ копирует repo_user_area/* (equ, exp, exe, sbr, fnc, ...) в корень
 ├─ make -f exe/Makefile clean
 └─ exe/as_exe -m flux_feqis -v aug34954 -s 4 -e 5      ← сразу первый запуск
       └─ as_exe.py → exe/Build
             ├─ source platform/env.<platform> + exe/astra_rc[_gcc|_ifx]  (пути к внешним библиотекам в $ASTRA_EXT)
             ├─ make exe/Makefile   → pyparse генерирует Fortran из модели, компилирует libastra/feqis/nbi/graph/... и линкует bin/<model>_gui.exe
             ├─ make exe/Makexpr    → вспомогательные IPC-процессы TGLF / QuaLiKiZ / NEO (нужен MPI)
             └─ запуск bin/<model>_gui.exe (X11-окно)
```

Ключевая идея: ASTRA **компилируется при каждом запуске модели**. Сторонние модули (json-fortran, SPIDER, RABBIT, TORBEAM, TGLF, NEO, QuaLiKiZ, ...) предполагаются **заранее собранными** в общем каталоге `$ASTRA_EXT` конкретной лаборатории. Их собирает `modules_install.sh` / `modules_inst_loc.sh`.

## Почему не ставится: найденные причины

1. **Нет платформы для обычного Linux.** На любом «чужом» хосте `get_platform` возвращает `mypc`, но файла `platform/env.mypc` нет. В итоге `source` падает, а `FC`, `MKL_LIBDIR` и `ASTRA_EXT` остаются пустыми. Это и есть «русские токамаки забыли добавить»: env-файлы есть только для конкретных кластеров (IPP, ITER, GA, MIT, ...).
2. **Поддержка обычного Linux у разработчиков недоделана.** `modules_inst_loc.sh` (добавлен 2026-08-10, «More portability for MacOS, Linux PC») ветвится по `platform == "linux"` (ставит пакеты через apt), однако `get_platform` никогда не возвращает `linux`, и `env.linux` в репозитории нет. Для macOS (`env.darwin` + brew + OpenBLAS) путь доведён до конца, для Ubuntu — нет.
3. **Жёсткие зависимости в линковке:**
   - **json-fortran** подключается безусловно (`JSON_LIB`/`JSON_INC` в `astra_rc*`), а в Ubuntu такого пакета нет: собирать из исходников (это делает `modules_inst_loc.sh`, публичный GitHub);
   - **MKL** подключается безусловно (`libmkl_*.a`, статически). Альтернатива OpenBLAS включается только при `PLATFORM == darwin`. На Ubuntu есть `libmkl-dev` (2020.4) или Intel oneAPI; второй вариант — патч `Makefile`/`Makexpr`, чтобы OpenBLAS выбирался и для Linux.
4. **Закрытые модули.** Без доступа к IPP недоступны RABBIT (NBI), TORBEAM (ECRH), SPIDER (свободная граница), STRAHL. С мая 2026 (коммит `2ec99235`) в коде есть заглушки на случай их отсутствия (`sbr/rabbit.F90`, `sbr/torba.F90`, `#ifdef SPIDER`), то есть **сборка без них возможна**. Но:
   - `install.sh` без флага `-safe` не отключает `RABBIT`/`TORBA` в тестовой модели `flux_feqis`, и первый автозапуск сразу идёт на модели, требующей закрытых кодов;
   - `Makexpr` по умолчанию собирает `tglfi`, `qlki`, `neo`. `QLK_LIB` прописан без проверки существования файла, поэтому без QuaLiKiZ/GACODE и MPI сборка IPC-части падает (`Build` завершается с кодом 2 ещё до запуска модели). Выход: `install.sh -safe` и ответить `n` на QLK/NEO/TGLF или вычистить цели в `Makexpr`.
5. **SSH-URL в скриптах** (`git@...`) даже для публичных репозиториев: без настроенных ключей клоны падают.
6. **Компилятор.** «Родной» для кода — Intel (`ifort`/`ifx`, флаги `-r8 -no-prec-div`). Поддержку gfortran довели в 2026 (`astra_rc_gcc`, `env.tok_gcc`, коммиты «gfortran working», «double precision everywhere»), но проверяли её только на кластере MPCDF с gcc 14. На Ubuntu 24.04 по умолчанию gfortran 13, совместимость не проверена.

Мелочь со скриншотов A7.02: `rm: No match` (csh-синтаксис в скриптах) и `cp: cannot create '.tsk/Intel/...exe'` (нет каталога `.tsk/Intel`) — на работу не влияют, но говорят о том же: скрипты писались под конкретную машину.

## Что нужно для минимальной рабочей сборки на Ubuntu

Минимальный состав: модели `showdata`/`readme` и `flux_feqis` с отключёнными RABBIT/TORBEAM.

| Компонент | Откуда | Обязателен |
|---|---|---|
| gfortran, gcc, make, cmake | apt | да |
| python3 + numpy + scipy | apt / venv | да (парсер моделей `pyparse`, `as_exe.py`) |
| libx11-dev + X-сервер (на Windows — WSLg) | apt | да для GUI; без него только batch (`-batch`) |
| MKL (`libmkl-dev`) **или** OpenBLAS + патч Makefile | apt | да |
| json-fortran 9.0.2 | сборка из GitHub | да |
| NetCDF (+fortran) | apt | только для RABBIT/вывода в CDF |
| OpenMPI | apt | только для TGLF/NEO/QuaLiKiZ |
| TGLF/NEO (GACODE), QuaLiKiZ | публичные, сборка | нет |
| RABBIT, TORBEAM, SPIDER, STRAHL | закрыты (IPP) | нет; без них нет NBI/ECRH/свободной границы |

Вывод: базовый запуск на Ubuntu **реалистичен без закрытых кодов**. Основная работа — написать `platform/env.mypc` (или `env.linux`) для gfortran + apt-пакетов и собрать json-fortran в локальный `$ASTRA_EXT`.

## План

1. Запросить у Никиты `install A8.txt`, файлы `readme`/`showdata` и эталонные выводы A7.02.
2. Поднять чистый Ubuntu 24.04 (контейнер/VM или WSL) и пройти установку с нуля, фиксируя каждый шаг.
3. Учесть ловушку с компилятором: `Build` берёт `exe/astra_rc` (пути Intel), если не задан `ASTRA_COMPILER`. А если задать `ASTRA_COMPILER=gcc` (или `as_exe -c gcc`), `get_platform` превращает имя платформы в `mypc_gcc`, и тогда нужен файл `env.mypc_gcc`. Написать env-файл для gfortran: `FC=gfortran`, `FC_FLAGS="-O0 -w -fPIC -fopenmp"`, `MKL_LIBDIR=/usr/lib/x86_64-linux-gnu`, `ASTRA_EXT=$AWD/astra_modules`. Привести пути в `astra_rc_gcc` к реальной раскладке `$ASTRA_EXT`: в `modules_inst_loc.sh` раскладка без версий, а `astra_rc_gcc` ждёт `json/gcc_9.0.2` и т. п.
4. Собрать json-fortran; выключить в `Makexpr` цели, для которых нет библиотек; запустить `showdata`/`readme`, затем `flux_feqis`/`aug34954` без RABBIT/TORBA и сравнить с регрессионными CDF из `regressions/`.
5. Оформить runbook «чистый Ubuntu → работающая ASTRA 8», при возможности одним скриптом и/или Dockerfile.
6. Отдельно: имеет смысл предложить мейнтейнерам (G. Tardini) PR с generic-платформой `linux`/`mypc`, раз у них уже есть такая заготовка.

## Открытые вопросы к лаборатории

- Что именно уже сделано и на чём остановились (`install A8.txt`)?
- Нужны ли на первом этапе NBI (RABBIT) и ECRH (TORBEAM)? Если да, понадобится доступ через сотрудничество с IPP.
- Целевая машина: WSL на Windows (как на скриншоте A7.02) или отдельный Linux-ПК? От этого зависит, как выводить X11.
- Есть ли лицензия/желание использовать Intel oneAPI (бесплатный) вместо gfortran?

## Итог (2026-10-01, вечер)

Сборка на Ubuntu 24.04 (WSL2) получилась без закрытых модулей. Установщик и документация лежат в отдельной папке `/home/egor/Work/astra-installer/`: `install_astra8.sh`, `README.md` (коротко), `docs/DETAILED.md` (подробно). Тестовый расчёт `flux_feqis`/`aug34954` доходит до `ASTRA normal exit`, GUI-окно работает через WSLg. Случай `fbe`/`AUG33040_2500` совпадает с эталоном IPP на уровне ~0.1 % по току и q.

Кроме перечисленного выше, при сборке нашлись ещё проблемы: CRLF при `core.autocrlf=true`, строки длиннее 132 символов (`-ffree-line-length-none`), `--as-needed` ломает линковку X11 (`-Wl,--no-as-needed`), парсер U-строк в формате A7 и переподписка потоков OpenBLAS внутри OpenMP. `showdata`/`readme` не запускается без U-файлов `00000.tes`/`00000.tis` из установки A7 лаборатории.
