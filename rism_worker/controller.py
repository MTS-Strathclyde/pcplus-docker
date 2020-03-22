import logging
from .rism import run_rism
from .models import CalculationRequest, CalculationResponse

logger = logging.getLogger(__name__)


def calculate(calc_request) -> CalculationResponse:
    results = None
    converged = False
    try:
        results = run_rism(calc_request)
        converged = True
    except Exception as e:
        logger.exception(e)
    return CalculationResponse(converged, results)
