# Документация

## Установка (для пользователей)

| документ | для кого |
|---|---|
| [`../installer/README.md`](../installer/README.md) | быстрая установка на Ubuntu / WSL |
| [`../installer/docker/README.md`](../installer/docker/README.md) | то же в Docker |
| [`../installer/docs/DETAILED.md`](../installer/docs/DETAILED.md) | как устроена сборка ASTRA, что не работает «из коробки» и как это исправлено |

## Изучение ASTRA

| документ | о чём |
|---|---|
| [`astra/platforms.md`](astra/platforms.md) | конфиги площадок `platform/env.*`: как ASTRA выбирает площадку, чем отличаются IPP, ITER, GA, NERSC и др. |

Сюда же складываем разборы: устройство моделей (`equ/`), парсер `pyparse/`, FEQIS, вывод и т.д.

## Форк и процесс

| документ | о чём |
|---|---|
| [`local-fork.md`](local-fork.md) | `astra-src/`: ветки, синхронизация с оригиналом, лицензия, журнал правок |

## Связанные проекты

- **tokomak** — `/home/egor/Work/tokomak` (из Windows: `\\wsl.localhost\Ubuntu\home\egor\Work\tokomak`),
  наш внутренний проект лаборатории (архив «Глобус-М2»). Там много контекста про установку и про то,
  как всё устроено в лаборатории: сервер, доступ, деплой. Пока не разобран, ссылка на будущее.

## Исследования

| документ | о чём |
|---|---|
| [`research/astra-a8-install-research-2026-10-01.md`](research/astra-a8-install-research-2026-10-01.md) | исходный запрос лаборатории и разбор, почему ASTRA 8 не ставится на обычную Ubuntu |
