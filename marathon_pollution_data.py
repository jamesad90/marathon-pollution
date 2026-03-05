"""
Get pollution data for World Marathon Majors
Retrieves PM2.5, PM10, and NO2 data for:
- Race day
- 7 days before race day

Uses Copernicus CAMS data
"""

import pandas as pd
import cdsapi
import xarray as xr
import numpy as np
from datetime import datetime, timedelta
import time
import json
import os

# Marathon city coordinates
MARATHON_CITIES = {
    'Boston Marathon': {
        'city': 'Boston',
        'country': 'USA',
        'latitude': 42.3601,
        'longitude': -71.0589
    },
    'London Marathon': {
        'city': 'London',
        'country': 'UK',
        'latitude': 51.5074,
        'longitude': -0.1278
    },
    'Berlin Marathon': {
        'city': 'Berlin',
        'country': 'Germany',
        'latitude': 52.5200,
        'longitude': 13.4050
    },
    'Chicago Marathon': {
        'city': 'Chicago',
        'country': 'USA',
        'latitude': 41.8781,
        'longitude': -87.6298
    },
    'New York City Marathon': {
        'city': 'New York',
        'country': 'USA',
        'latitude': 40.7128,
        'longitude': -74.0060
    },
    'Tokyo Marathon': {
        'city': 'Tokyo',
        'country': 'Japan',
        'latitude': 35.6762,
        'longitude': 139.6503
    }
}

# Marathon typical race months (for determining dates)
# We'll need to look these up, but typical patterns:
MARATHON_MONTHS = {
    'Boston Marathon': 4,  # April (3rd Monday)
    'London Marathon': 4,  # April
    'Berlin Marathon': 9,  # September
    'Chicago Marathon': 10,  # October
    'New York City Marathon': 11,  # November (1st Sunday)
    'Tokyo Marathon': 3  # March (was February pre-2020)
}


def get_marathon_dates():
    """
    Get actual marathon race dates for each year
    Returns DataFrame with columns: year, marathon, date
    """
    dates_data = []

    # Known marathon dates (this would ideally be more complete)
    # For demonstration, I'll add known dates from our data
    marathon_dates = {
        2010: {
            'Boston Marathon': '2010-04-19',
            'London Marathon': '2010-04-25',
            'Berlin Marathon': '2010-09-26',
            'Chicago Marathon': '2010-10-10',
            'New York City Marathon': '2010-11-07',
            'Tokyo Marathon': '2010-02-28'
        },
        2011: {
            'Boston Marathon': '2011-04-18',
            'London Marathon': '2011-04-17',
            'Berlin Marathon': '2011-09-25',
            'Chicago Marathon': '2011-10-09',
            'New York City Marathon': '2011-11-06',
            'Tokyo Marathon': '2011-02-27'
        },
        2012: {
            'Boston Marathon': '2012-04-16',
            'London Marathon': '2012-04-22',
            'Berlin Marathon': '2012-09-30',
            'Chicago Marathon': '2012-10-07',
            'New York City Marathon': None,  # Cancelled
            'Tokyo Marathon': '2012-02-26'
        },
        2013: {
            'Boston Marathon': '2013-04-15',
            'London Marathon': '2013-04-21',
            'Berlin Marathon': '2013-09-29',
            'Chicago Marathon': '2013-10-13',
            'New York City Marathon': '2013-11-03',
            'Tokyo Marathon': '2013-02-24'
        },
        2014: {
            'Boston Marathon': '2014-04-21',
            'London Marathon': '2014-04-13',
            'Berlin Marathon': '2014-09-28',
            'Chicago Marathon': '2014-10-12',
            'New York City Marathon': '2014-11-02',
            'Tokyo Marathon': '2014-02-23'
        },
        2015: {
            'Boston Marathon': '2015-04-20',
            'London Marathon': '2015-04-26',
            'Berlin Marathon': '2015-09-27',
            'Chicago Marathon': '2015-10-11',
            'New York City Marathon': '2015-11-01',
            'Tokyo Marathon': '2015-02-22'
        },
        2016: {
            'Boston Marathon': '2016-04-18',
            'London Marathon': '2016-04-24',
            'Berlin Marathon': '2016-09-25',
            'Chicago Marathon': '2016-10-09',
            'New York City Marathon': '2016-11-06',
            'Tokyo Marathon': '2016-02-28'
        },
        2017: {
            'Boston Marathon': '2017-04-17',
            'London Marathon': '2017-04-23',
            'Berlin Marathon': '2017-09-24',
            'Chicago Marathon': '2017-10-08',
            'New York City Marathon': '2017-11-05',
            'Tokyo Marathon': '2017-02-26'
        },
        2018: {
            'Boston Marathon': '2018-04-16',
            'London Marathon': '2018-04-22',
            'Berlin Marathon': '2018-09-16',
            'Chicago Marathon': '2018-10-07',
            'New York City Marathon': '2018-11-04',
            'Tokyo Marathon': '2018-02-25'
        },
        2019: {
            'Boston Marathon': '2019-04-15',
            'London Marathon': '2019-04-28',
            'Berlin Marathon': '2019-09-29',
            'Chicago Marathon': '2019-10-13',
            'New York City Marathon': '2019-11-03',
            'Tokyo Marathon': '2019-03-03'
        },
        2020: {
            'Boston Marathon': None,  # Cancelled
            'London Marathon': None,  # Cancelled/Virtual
            'Berlin Marathon': None,  # Cancelled
            'Chicago Marathon': None,  # Cancelled
            'New York City Marathon': None,  # Cancelled
            'Tokyo Marathon': '2020-03-01'  # Elite only
        },
        2021: {
            'Boston Marathon': '2021-10-11',
            'London Marathon': '2021-10-03',
            'Berlin Marathon': '2021-09-26',
            'Chicago Marathon': '2021-10-10',
            'New York City Marathon': '2021-11-07',
            'Tokyo Marathon': '2021-03-07'
        },
        2022: {
            'Boston Marathon': '2022-04-18',
            'London Marathon': '2022-10-02',
            'Berlin Marathon': '2022-09-25',
            'Chicago Marathon': '2022-10-09',
            'New York City Marathon': '2022-11-06',
            'Tokyo Marathon': '2022-03-06'
        },
        2023: {
            'Boston Marathon': '2023-04-17',
            'London Marathon': '2023-04-23',
            'Berlin Marathon': '2023-09-24',
            'Chicago Marathon': '2023-10-08',
            'New York City Marathon': '2023-11-05',
            'Tokyo Marathon': '2023-03-05'
        },
        2024: {
            'Boston Marathon': '2024-04-15',
            'London Marathon': '2024-04-21',
            'Berlin Marathon': '2024-09-29',
            'Chicago Marathon': '2024-10-13',
            'New York City Marathon': '2024-11-03',
            'Tokyo Marathon': '2024-03-03'
        },
        2025: {
            'Boston Marathon': '2025-04-21',
            'London Marathon': '2025-04-27',
            'Berlin Marathon': '2025-09-28',
            'Chicago Marathon': '2025-10-12',
            'New York City Marathon': '2025-11-02',
            'Tokyo Marathon': '2025-03-02'
        }
    }

    for year, marathons in marathon_dates.items():
        for marathon_name, date_str in marathons.items():
            if date_str:  # Skip cancelled races
                dates_data.append({
                    'year': year,
                    'marathon': marathon_name,
                    'race_date': date_str,
                    'city': MARATHON_CITIES[marathon_name]['city'],
                    'country': MARATHON_CITIES[marathon_name]['country'],
                    'latitude': MARATHON_CITIES[marathon_name]['latitude'],
                    'longitude': MARATHON_CITIES[marathon_name]['longitude']
                })

    return pd.DataFrame(dates_data)


def download_daily_pollution_data(client, year, month, day, bbox, output_file):
    """
    Download daily pollution data from CAMS

    Parameters:
    - client: CDS API client
    - year, month, day: Date components
    - bbox: [north, west, south, east] bounding box
    - output_file: path to save the netCDF file
    """

    if year < 2003:
        print(f"Year {year} is before CAMS data availability (2003)")
        return None

    request = {
        'format': 'netcdf',
        'variable': [
            'particulate_matter_2.5um',
            'particulate_matter_10um',
            'nitrogen_dioxide'
        ],
        'year': str(year),
        'month': str(month).zfill(2),
        'day': str(day).zfill(2),
        'time': ['00:00', '06:00', '12:00', '18:00'],  # 4 times per day
        'area': bbox,
        'product_type': 'reanalysis'
    }

    try:
        print(f"Downloading CAMS data for {year}-{month:02d}-{day:02d}...")
        client.retrieve(
            'cams-global-reanalysis-eac4',
            request,
            output_file
        )
        print(f"Data saved to {output_file}")
        return output_file
    except Exception as e:
        print(f"Error downloading data for {year}-{month:02d}-{day:02d}: {e}")
        return None


def extract_city_daily_data(nc_file, latitude, longitude, date_str):
    """
    Extract pollution data for a specific city from daily netCDF file

    Returns dictionary with average values for the day
    """

    try:
        ds = xr.open_dataset(nc_file)

        # Find nearest grid point
        city_point = ds.sel(
            latitude=latitude,
            longitude=longitude,
            method='nearest'
        )

        # Calculate daily averages across all time points
        result = {
            'date': date_str,
            'latitude': latitude,
            'longitude': longitude
        }

        # Map variable names to our desired names
        var_mapping = {
            'pm2p5': 'PM2.5',
            'pm10': 'PM10',
            'no2': 'NO2'
        }

        for var_key, var_name in var_mapping.items():
            if var_key in ds.data_vars:
                daily_mean = float(city_point[var_key].mean())
                result[var_name] = daily_mean

        ds.close()
        return result

    except Exception as e:
        print(f"Error extracting data from {nc_file}: {e}")
        return None


def get_pollution_for_marathon(client, marathon_row, output_dir='pollution_data'):
    """
    Get pollution data for a specific marathon (race day + 7 days before)

    Parameters:
    - marathon_row: Row from marathon dates DataFrame
    - output_dir: Directory to save netCDF files

    Returns:
    - DataFrame with pollution data for race day and 7 days before
    """

    os.makedirs(output_dir, exist_ok=True)

    race_date = datetime.strptime(marathon_row['race_date'], '%Y-%m-%d')

    # Calculate 7 days before
    week_before = race_date - timedelta(days=7)

    # Get bounding box for the city (small area around the city)
    lat = marathon_row['latitude']
    lon = marathon_row['longitude']
    bbox = [lat + 2, lon - 2, lat - 2, lon + 2]  # [N, W, S, E]

    results = []

    # Get data for both dates
    for target_date in [week_before, race_date]:
        date_type = 'race_day' if target_date == race_date else '7_days_before'

        nc_file = os.path.join(
            output_dir,
            f"{marathon_row['marathon'].replace(' ', '_')}_{target_date.strftime('%Y%m%d')}.nc"
        )

        # Download data
        downloaded_file = download_daily_pollution_data(
            client,
            target_date.year,
            target_date.month,
            target_date.day,
            bbox,
            nc_file
        )

        if downloaded_file:
            # Extract city data
            city_data = extract_city_daily_data(
                downloaded_file,
                lat,
                lon,
                target_date.strftime('%Y-%m-%d')
            )

            if city_data:
                city_data['year'] = marathon_row['year']
                city_data['marathon'] = marathon_row['marathon']
                city_data['city'] = marathon_row['city']
                city_data['country'] = marathon_row['country']
                city_data['date_type'] = date_type
                city_data['race_date'] = marathon_row['race_date']
                results.append(city_data)

        # Be polite to the API
        time.sleep(2)

    return pd.DataFrame(results) if results else pd.DataFrame()


def process_all_marathons(start_year=2010, end_year=2025):
    """
    Process all World Marathon Majors for specified years
    """

    print("=" * 80)
    print("World Marathon Majors Pollution Data Collection")
    print("=" * 80)
    print()

    # Get marathon dates
    marathon_dates = get_marathon_dates()
    marathon_dates = marathon_dates[
        (marathon_dates['year'] >= start_year) &
        (marathon_dates['year'] <= end_year)
    ]

    print(f"Processing {len(marathon_dates)} marathon races from {start_year} to {end_year}")
    print()

    # Setup CDS API
    print("Setting up Copernicus CDS API...")
    try:
        client = cdsapi.Client()
    except Exception as e:
        print(f"Error setting up CDS API: {e}")
        print("Make sure you have ~/.cdsapirc configured with your API key")
        return None

    # Process each marathon
    all_results = []

    for idx, row in marathon_dates.iterrows():
        print(f"\n[{idx+1}/{len(marathon_dates)}] Processing {row['marathon']} {row['year']}...")
        print(f"Race date: {row['race_date']}")

        try:
            marathon_pollution = get_pollution_for_marathon(client, row)
            if not marathon_pollution.empty:
                all_results.append(marathon_pollution)
                print(f"✓ Successfully collected pollution data")
            else:
                print(f"✗ No pollution data collected")
        except Exception as e:
            print(f"✗ Error processing marathon: {e}")

    # Combine all results
    if all_results:
        results_df = pd.concat(all_results, ignore_index=True)

        # Save to CSV
        output_file = f'marathon_pollution_{start_year}_{end_year}.csv'
        results_df.to_csv(output_file, index=False)

        print("\n" + "=" * 80)
        print(f"✓ Results saved to {output_file}")
        print(f"✓ Total records: {len(results_df)}")
        print("=" * 80)

        # Display summary
        print("\nSummary by date type:")
        print(results_df.groupby('date_type').size())

        print("\nSample pollution data:")
        print(results_df[['marathon', 'year', 'date_type', 'PM2.5', 'PM10', 'NO2']].head(10))

        return results_df
    else:
        print("\n✗ No data was successfully retrieved")
        return pd.DataFrame()


if __name__ == "__main__":
    # Test with a small sample first
    print("Starting marathon pollution data collection...")
    print("This will take some time due to API rate limits")
    print()

    # Start with just 2023-2024 for testing
    results = process_all_marathons(start_year=2023, end_year=2024)

    # Uncomment to process all years:
    # results = process_all_marathons(start_year=2010, end_year=2025)
