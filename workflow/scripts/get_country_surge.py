from pathlib import Path
import sys
import logging, traceback

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

global_haz = snakemake.input.global_haz
climate_scenario = snakemake.wildcards.climate_scenario
ssp = snakemake.wildcards.ssp
slr_year = snakemake.wildcards.slr_year
cc = coco.CountryConverter()
countries = snakemake.params.countries

countries_iso_num = cc.convert(countries, to="ISOnumeric")

logger.info(f"Splitting Surge events in countries for {climate_scenario}/{ssp}/{slr_year}slr from {global_haz}")
haz = Hazard.from_hdf5(global_haz)
haz.centroid.set_region_id()
for cnt_id, country in zip(countries_iso_num, countries):
    haz_cnt = haz.select(reg_id=cnt_id)
    filename = f'surge/{climate_scenario}/surge_{country}_{ssp}_{slr_year}slr_{climate_scenario}.hdf5'
    if haz_cnt:
        haz_cnt.write_hdf5(filename)
    else:
        Path(filename).touch()
