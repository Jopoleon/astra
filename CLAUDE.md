# CLAUDE.md

Репозиторий `Jopoleon/astra` — наша обвязка вокруг транспортного кода ASTRA 8 (документация,
установщик для Ubuntu, будущий apt-пакет). Язык документации — русский.

## Главное

- `astra-src/` — локальный форк исходников ASTRA (оригинал MPCDF, LGPL v2.1+). Папка в `.gitignore`,
  у неё свой git: remote `upstream` = MPCDF, `main` = копия оригинала (не трогать),
  `local` = наша рабочая ветка. Подробно: `docs/local-fork.md`.
- Правки исходников ASTRA делаем только в `astra-src/` на ветке `local`, коммитим там же.
  Значимые правки и удаления отмечаем в журнале в `docs/local-fork.md`.
- Не пушить `astra-src/` в GitHub: в истории есть файлы >100 МБ.
- Соблюдать LGPL: не удалять `LICENSE` и копирайты оригинала.

## Связанное

- `/home/egor/Work/astra-installer/` — текущий установщик ASTRA 8 для Ubuntu (`install_astra8.sh`,
  `docs/DETAILED.md` — как устроена сборка ASTRA и какие патчи нужны для Ubuntu).
- `/home/egor/Work/tokomak/` — основной проект лаборатории (архив «Глобус-М2»); ASTRA там описана
  как побочная задача в `docs/context/project-context.md`.
