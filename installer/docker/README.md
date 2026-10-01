# ASTRA 8 в Docker

Тот же `install_astra8.sh`, но внутри контейнера Ubuntu 24.04: в систему ставить ничего не нужно, кроме Docker. Сборка образа сама прогоняет тестовый расчёт, поэтому если образ собрался, ASTRA в нём работает.

## Что нужно

- Docker: на Windows — Docker Desktop с включённой интеграцией WSL (Settings → Resources → WSL integration → ваша Ubuntu); на Linux — `docker.io` или Docker Engine;
- ~3 ГБ на диске под образ;
- для окна с графиками: Windows 11 + WSL2 (WSLg) или рабочий стол Linux.

## Запуск

Из папки `installer/` нашего репозитория:

```bash
# первый запуск сам соберёт образ "astra8" (5-10 минут, один раз)
docker/astra-docker.sh -m flux_feqis -v aug34954 -s 4 -e 5          # с окном графиков
docker/astra-docker.sh -m flux_feqis -v aug34954 -s 4 -e 5 -batch   # без окна
docker/astra-docker.sh                                              # повторить последний расчёт
docker/astra-docker.sh shell                                        # bash внутри контейнера
```

Аргументы такие же, как у `astra8` (см. [`../README.md`](../README.md)).

## Где файлы

Папка `~/astra-data` на компьютере подключается в контейнер как `/data`:

| на компьютере | что там |
|---|---|
| `~/astra-data/equ/` | модели — сюда класть свои |
| `~/astra-data/exp/` | файлы эксперимента |
| `~/astra-data/udb/` | U-файлы |
| `~/astra-data/ncdf_out/` | результаты (`<данные><модель>.CDF`) |

При первом запуске туда копируются примеры из образа. Ваши файлы не перезаписываются, контейнер после расчёта удаляется, а данные остаются. Другая папка: `ASTRA_DATA=/путь docker/astra-docker.sh ...`.

## Сборка вручную

```bash
docker build -f docker/Dockerfile -t astra8 \
  --build-arg USER_UID=$(id -u) --build-arg USER_GID=$(id -g) .
```

Параметры сборки (`--build-arg`):

| параметр | по умолчанию | что это |
|---|---|---|
| `UBUNTU_VERSION` | `24.04` | базовый образ (`22.04` тоже проверен) |
| `ASTRA_REF` | `e66f6ba1...` (проверенный) | версия ASTRA: полный хэш коммита, `main` или тег |
| `BLAS` | `openblas` | `openblas` или `mkl` |
| `USER_UID`, `USER_GID` | `1000` | чтобы файлы в `~/astra-data` принадлежали вам (`astra-docker.sh` подставляет сам) |

Обновить ASTRA до свежей версии: `docker build ... --build-arg ASTRA_REF=main -t astra8 .`

## Файлы

| файл | назначение |
|---|---|
| `Dockerfile` | образ: Ubuntu + пользователь `astra` + `install_astra8.sh` |
| `entrypoint.sh` | внутри контейнера: подключает `/data`, запускает `astra8` |
| `astra-docker.sh` | на компьютере: собирает образ при необходимости, пробрасывает окно X11 и папку данных |
