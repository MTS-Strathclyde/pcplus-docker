import flask
import logging

app = flask.Flask(__name__)
app.config.from_object('config.Config')

if __name__ != '__main__':
    gunicorn_logger = logging.getLogger('gunicorn.error')
    app.logger.handlers = gunicorn_logger.handlers
    app.logger.setLevel(gunicorn_logger.level)

from . import routes  # noqa
