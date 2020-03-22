#!/bin/sh

# Putting gunicorn command directly into docker resulted in weird errrors
# Once I wrapped it with a shell script it started to work fine.

exec gunicorn --workers 3 -k gevent --bind 0.0.0.0:3002 wsgi:app \
    --max-requests 0 --timeout 3600 --keep-alive 3600 --log-level ${LOG_LEVEL:=info}
