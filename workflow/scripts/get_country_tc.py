import sys
import logging, traceback

import numpy as np
import country_converter as coco

from climada.hazard import Hazard

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())

def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(
        "".join(
            [
                "Uncaught exception: ",
                *traceback.format_exception(exc_type, exc_value, exc_traceback),
            ]
        )
    )



# Install exception handler
sys.excepthook = handle_exception

climate_scenario = snakemake.wildcards.climate_scenario
global_haz = snakemake.input.global_haz
cc = coco.CountryConverter()
countries = snakemake.params.countries

countries_iso = cc.convert(countries, to="ISOnumeric")
logger.info(f"Splitting TC events in countries for {climate_scenario} from {global_haz}")
haz = Hazard.from_hdf5(global_haz)
for country, cnt_id in zip(countries_iso, countries):
    haz_cnt = haz.select(reg_id=cnt_id)
    filename = f'tropcyc/{climate_scenario}/tropcyc_{country}_{climate_scenario}.hdf5'
    haz_cnt.write_hdf5(filename)