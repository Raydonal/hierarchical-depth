"""Fit the mixed-Poisson depth model to the OpenFlights route network.

Data: `data/openflights-airports.dat` and `data/openflights-routes.dat`,
snapshots of the OpenFlights Airport and Route Databases
(https://github.com/jpatokal/openflights), Open Database License (ODbL) --
free to use with attribution, no registration, no notify-on-publication
requirement (checked 2026-09-05; this is the reason it was picked over the
CAIDA AS-relationships dataset, which requires accepting an Acceptable Use
Agreement).

Unlike `data/email-Eu-core*` (a single flat 42-department partition),
airports carry a genuine THREE-level geographic nesting: continent > country
> airport. Continent is not in the raw data -- it's added here via a
hand-built country-name -> continent lookup (`COUNTRY_TO_CONTINENT` below),
covering every country string that actually appears in this snapshot. This
lookup is this script's own addition, not part of the original dataset.

IMPORTANT CAVEAT (the same one that applies to email-Eu-core): the estimated depth N from `fit_homogeneous_network`
is an *effective statistical depth* of the degree-generating mechanism, not a
literal count of geographic levels. Continent/country nesting is a plausible
*source* of overdispersion in an airport's route count (a hub in a well-
connected continent/country tends to have more routes), but finding
N_round != 1 here would not by itself prove the geographic hierarchy is
*mechanistically* what the cascade model describes -- it would only show the
degree sequence is statistically consistent with *some* depth-N mixture.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

SEED = 20260905
SYSTEM = "OpenFlights"

# Every country string appearing in data/openflights-airports.dat (237 as of
# the 2026-09-05 snapshot), mapped to one of six continents (Antarctica
# folded into "Oceania-adjacent" territories are kept literal). Built by hand
# for this script -- not sourced from a package -- so double-check before
# reuse elsewhere.
COUNTRY_TO_CONTINENT = {
    "Afghanistan": "Asia", "Albania": "Europe", "Algeria": "Africa",
    "American Samoa": "Oceania", "Angola": "Africa", "Anguilla": "North America",
    "Antarctica": "Antarctica", "Antigua and Barbuda": "North America",
    "Argentina": "South America", "Armenia": "Asia", "Aruba": "North America",
    "Australia": "Oceania", "Austria": "Europe", "Azerbaijan": "Asia",
    "Bahamas": "North America", "Bahrain": "Asia", "Bangladesh": "Asia",
    "Barbados": "North America", "Belarus": "Europe", "Belgium": "Europe",
    "Belize": "North America", "Benin": "Africa", "Bermuda": "North America",
    "Bhutan": "Asia", "Bolivia": "South America",
    "Bosnia and Herzegovina": "Europe", "Botswana": "Africa", "Brazil": "South America",
    "British Indian Ocean Territory": "Asia", "British Virgin Islands": "North America",
    "Brunei": "Asia", "Bulgaria": "Europe", "Burkina Faso": "Africa",
    "Burma": "Asia", "Burundi": "Africa", "Cambodia": "Asia", "Cameroon": "Africa",
    "Canada": "North America", "Cape Verde": "Africa", "Cayman Islands": "North America",
    "Central African Republic": "Africa", "Chad": "Africa", "Chile": "South America",
    "China": "Asia", "Christmas Island": "Oceania", "Cocos (Keeling) Islands": "Oceania",
    "Colombia": "South America", "Comoros": "Africa", "Congo (Brazzaville)": "Africa",
    "Congo (Kinshasa)": "Africa", "Cook Islands": "Oceania", "Costa Rica": "North America",
    "Cote d'Ivoire": "Africa", "Croatia": "Europe", "Cuba": "North America",
    "Cyprus": "Europe", "Czech Republic": "Europe", "Denmark": "Europe",
    "Djibouti": "Africa", "Dominica": "North America", "Dominican Republic": "North America",
    "East Timor": "Asia", "Ecuador": "South America", "Egypt": "Africa",
    "El Salvador": "North America", "Equatorial Guinea": "Africa", "Eritrea": "Africa",
    "Estonia": "Europe", "Ethiopia": "Africa", "Falkland Islands": "South America",
    "Faroe Islands": "Europe", "Fiji": "Oceania", "Finland": "Europe", "France": "Europe",
    "French Guiana": "South America", "French Polynesia": "Oceania", "Gabon": "Africa",
    "Gambia": "Africa", "Georgia": "Asia", "Germany": "Europe", "Ghana": "Africa",
    "Gibraltar": "Europe", "Greece": "Europe", "Greenland": "North America",
    "Grenada": "North America", "Guadeloupe": "North America", "Guam": "Oceania",
    "Guatemala": "North America", "Guernsey": "Europe", "Guinea": "Africa",
    "Guinea-Bissau": "Africa", "Guyana": "South America", "Haiti": "North America",
    "Honduras": "North America", "Hong Kong": "Asia", "Hungary": "Europe",
    "Iceland": "Europe", "India": "Asia", "Indonesia": "Asia", "Iran": "Asia",
    "Iraq": "Asia", "Ireland": "Europe", "Isle of Man": "Europe", "Israel": "Asia",
    "Italy": "Europe", "Jamaica": "North America", "Japan": "Asia", "Jersey": "Europe",
    "Johnston Atoll": "Oceania", "Jordan": "Asia", "Kazakhstan": "Asia", "Kenya": "Africa",
    "Kiribati": "Oceania", "Kuwait": "Asia", "Kyrgyzstan": "Asia", "Laos": "Asia",
    "Latvia": "Europe", "Lebanon": "Asia", "Lesotho": "Africa", "Liberia": "Africa",
    "Libya": "Africa", "Lithuania": "Europe", "Luxembourg": "Europe", "Macau": "Asia",
    "Macedonia": "Europe", "Madagascar": "Africa", "Malawi": "Africa", "Malaysia": "Asia",
    "Maldives": "Asia", "Mali": "Africa", "Malta": "Europe", "Marshall Islands": "Oceania",
    "Martinique": "North America", "Mauritania": "Africa", "Mauritius": "Africa",
    "Mayotte": "Africa", "Mexico": "North America", "Micronesia": "Oceania",
    "Midway Islands": "Oceania", "Moldova": "Europe", "Mongolia": "Asia",
    "Montenegro": "Europe", "Montserrat": "North America", "Morocco": "Africa",
    "Mozambique": "Africa", "Myanmar": "Asia", "Namibia": "Africa", "Nauru": "Oceania",
    "Nepal": "Asia", "Netherlands": "Europe", "Netherlands Antilles": "North America",
    "New Caledonia": "Oceania", "New Zealand": "Oceania", "Nicaragua": "North America",
    "Niger": "Africa", "Nigeria": "Africa", "Niue": "Oceania", "Norfolk Island": "Oceania",
    "North Korea": "Asia", "Northern Mariana Islands": "Oceania", "Norway": "Europe",
    "Oman": "Asia", "Pakistan": "Asia", "Palau": "Oceania", "Palestine": "Asia",
    "Panama": "North America", "Papua New Guinea": "Oceania", "Paraguay": "South America",
    "Peru": "South America", "Philippines": "Asia", "Poland": "Europe",
    "Portugal": "Europe", "Puerto Rico": "North America", "Qatar": "Asia",
    "Reunion": "Africa", "Romania": "Europe", "Russia": "Europe", "Rwanda": "Africa",
    "Saint Helena": "Africa", "Saint Kitts and Nevis": "North America",
    "Saint Lucia": "North America", "Saint Pierre and Miquelon": "North America",
    "Saint Vincent and the Grenadines": "North America", "Samoa": "Oceania",
    "Sao Tome and Principe": "Africa", "Saudi Arabia": "Asia", "Senegal": "Africa",
    "Serbia": "Europe", "Seychelles": "Africa", "Sierra Leone": "Africa",
    "Singapore": "Asia", "Slovakia": "Europe", "Slovenia": "Europe",
    "Solomon Islands": "Oceania", "Somalia": "Africa", "South Africa": "Africa",
    "South Korea": "Asia", "South Sudan": "Africa", "Spain": "Europe",
    "Sri Lanka": "Asia", "Sudan": "Africa", "Suriname": "South America",
    "Svalbard": "Europe", "Swaziland": "Africa", "Sweden": "Europe",
    "Switzerland": "Europe", "Syria": "Asia", "Taiwan": "Asia", "Tajikistan": "Asia",
    "Tanzania": "Africa", "Thailand": "Asia", "Togo": "Africa", "Tonga": "Oceania",
    "Trinidad and Tobago": "North America", "Tunisia": "Africa", "Turkey": "Asia",
    "Turkmenistan": "Asia", "Turks and Caicos Islands": "North America",
    "Tuvalu": "Oceania", "Uganda": "Africa", "Ukraine": "Europe",
    "United Arab Emirates": "Asia", "United Kingdom": "Europe",
    "United States": "North America", "Uruguay": "South America",
    "Uzbekistan": "Asia", "Vanuatu": "Oceania", "Venezuela": "South America",
    "Vietnam": "Asia", "Virgin Islands": "North America", "Wake Island": "Oceania",
    "Wallis and Futuna": "Oceania", "West Bank": "Asia", "Western Sahara": "Africa",
    "Yemen": "Asia", "Zambia": "Africa", "Zimbabwe": "Africa",
}


def main():
    from scripts import common, data_loaders

    rng_test = np.random.default_rng(SEED)
    rng_ci = np.random.default_rng([SEED, 1])
    airports = data_loaders.openflights_airports()
    ids, degrees = data_loaders.openflights_degrees()
    countries = [airports[i][0] if i in airports else "?" for i in ids]
    unmapped = set(countries) - set(COUNTRY_TO_CONTINENT) - {"?"}
    if unmapped:
        print(f"WARNING: {len(unmapped)} country strings with no continent mapping: {unmapped}")
    continents = {COUNTRY_TO_CONTINENT.get(c, "Unknown") for c in countries}
    nesting = {"n_airports": len(ids), "n_countries": len(set(countries) - {"?"}),
               "n_continents": len(continents - {"Unknown"}),
               "n_route_rows": int(degrees["out-degree"].sum())}
    print(f"OpenFlights route network: {nesting}")

    series = []
    for name, D in degrees.items():
        n_zero = int(np.sum(D == 0))
        series.append(common.analyze_degree_series(
            SYSTEM, name, D[D > 0], rng_test, rng_ci,
            n_zero=n_zero, zero_policy="dropped"))
    common.write_json("analyze_openflights", {
        "metadata": common.metadata(
            "scripts/analyze_openflights.py", SEED,
            rng_streams={"depth_test": "default_rng(SEED)",
                         "bootstrap_ci": "default_rng([SEED, 1])"},
            n_boot_test=common.N_BOOT_TEST, n_boot_ci=common.N_BOOT_CI,
            note=common.SELECTION_NOTE),
        "nesting": nesting,
        "series": series,
    })


if __name__ == "__main__":
    main()
