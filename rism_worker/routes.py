from typing import List

from flask import request, Response, json, jsonify, abort

from . import app
from .controller import calculate
from .models import CalculationRequest

import logging
logger = logging.getLogger(__name__)


@app.route("/health", methods=["GET"])
def health():
    return Response("OK", 200)


@app.route("/calculate", methods=["POST"])
def predict_route():
    try:
        content = request.get_json()
        logger.debug(content)
        calc_request = CalculationRequest(**content)
    except TypeError:
        logger.exception("Coulnd't parse request!")
        return jsonify({"error": "Couldn't parse the request, aborting!"}), 400
    try:
        response = calculate(calc_request)
        return jsonify(response.__dict__)
    except ValueError:
        logger.exception("Failed to run a job")
        return abort(500)
