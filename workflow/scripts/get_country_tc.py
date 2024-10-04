from pathlib import Path
import sys
import logging, traceback

import country_converter as coco

from climada.hazard import TropCyclone

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
tc_res_str = snakemake.params.tc_res_str
nsynth = snakemake.params.nsynth
start_period = snakemake.params.start_period
end_period = snakemake.params.end_period
countries = snakemake.params.countries

cc = coco.CountryConverter()

if countries != "all":
    countries_iso_num = cc.convert(countries, to="ISOnumeric")

logger.info(f"Splitting TC events in countries for {climate_scenario} from {global_haz}")
haz = TropCyclone.from_hdf5(global_haz)
haz.centroids.set_region_id()
if countries == "all":
    countries_tmp = cc.convert(haz.centroids.gdf["region_id"].unique(), to="ISO3", src="ISOnumeric")
    countries = [c for c in countries_tmp if c != "not found"]
    countries_iso_num = cc.convert(countries, to="ISOnumeric")
    
for cnt_id, country in zip(countries_iso_num, countries):
    haz_cnt = haz.select(reg_id=cnt_id)
    filename = f'tropcyc/{climate_scenario}/by_country/tropcyc_{tc_res_str}arcsec_{nsynth}synth_{country}_{start_period}_to_{end_period}_{climate_scenario}.hdf5'
    if haz_cnt:
        haz_cnt.write_hdf5(filename)
    else:
        Path(filename).touch()
