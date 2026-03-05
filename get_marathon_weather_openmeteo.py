"""
Get weather data for World Marathon Majors race days using Open-Meteo API

Open-Meteo provides free historical weather data (no API key needed)
Collects temperature, humidity, wind speed, and precipitation for race days
"""

import pandas as pd
import requests
from datetime import datetime
import time
import sys

# Marathon locations
MARATHON_LOCATIONS = {
    'Boston Marathon': {'city': 'Boston', 'country': 'USA', 'lat': 42.3601, 'lon': -71.0589},
    'London Marathon': {'city': 'London', 'country': 'UK', 'lat': 51.5074, 'lon': -0.1278},
    'Berlin Marathon': {'city': 'Berlin', 'country': 'Germany', 'lat': 52.5200, 'lon': 13.4050},
    'Chicago Marathon': {'city': 'Chicago', 'country': 'USA', 'lat': 41.8781, 'lon': -87.6298},
    'New York City Marathon': {'city': 'New York', 'country': 'USA', 'lat': 40.7128, 'lon': -74.0060},
    'Tokyo Marathon': {'city': 'Tokyo', 'country': 'Japan', 'lat': 35.6762, 'lon': 139.6503}
}

# Marathon race dates
MARATHON_DATES = {
    2010: {'Boston Marathon': '2010-04-19', 'London Marathon': '2010-04-25', 'Berlin Marathon': '2010-09-26', 'Chicago Marathon': '2010-10-10', 'New York City Marathon': '2010-11-07', 'Tokyo Marathon': '2010-02-28'},
    2011: {'Boston Marathon': '2011-04-18', 'London Marathon': '2011-04-17', 'Berlin Marathon': '2011-09-25', 'Chicago Marathon': '2011-10-09', 'New York City Marathon': '2011-11-06', 'Tokyo Marathon': '2011-02-27'},
    2012: {'Boston Marathon': '2012-04-16', 'London Marathon': '2012-04-22', 'Berlin Marathon': '2012-09-30', 'Chicago Marathon': '2012-10-07', 'Tokyo Marathon': '2012-02-26'},
    2013: {'Boston Marathon': '2013-04-15', 'London Marathon': '2013-04-21', 'Berlin Marathon': '2013-09-29', 'Chicago Marathon': '2013-10-13', 'New York City Marathon': '2013-11-03', 'Tokyo Marathon': '2013-02-24'},
    2014: {'Boston Marathon': '2014-04-21', 'London Marathon': '2014-04-13', 'Berlin Marathon': '2014-09-28', 'Chicago Marathon': '2014-10-12', 'New York City Marathon': '2014-11-02', 'Tokyo Marathon': '2014-02-23'},
    2015: {'Boston Marathon': '2015-04-20', 'London Marathon': '2015-04-26', 'Berlin Marathon': '2015-09-27', 'Chicago Marathon': '2015-10-11', 'New York City Marathon': '2015-11-01', 'Tokyo Marathon': '2015-02-22'},
    2016: {'Boston Marathon': '2016-04-18', 'London Marathon': '2016-04-24', 'Berlin Marathon': '2016-09-25', 'Chicago Marathon': '2016-10-09', 'New York City Marathon': '2016-11-06', 'Tokyo Marathon': '2016-02-28'},
    2017: {'Boston Marathon': '2017-04-17', 'London Marathon': '2017-04-23', 'Berlin Marathon': '2017-09-24', 'Chicago Marathon': '2017-10-08', 'New York City Marathon': '2017-11-05', 'Tokyo Marathon': '2017-02-26'},
    2018: {'Boston Marathon': '2018-04-16', 'London Marathon': '2018-04-22', 'Berlin Marathon': '2018-09-16', 'Chicago Marathon': '2018-10-07', 'New York City Marathon': '2018-11-04', 'Tokyo Marathon': '2018-02-25'},
    2019: {'Boston Marathon': '2019-04-15', 'London Marathon': '2019-04-28', 'Berlin Marathon': '2019-09-29', 'Chicago Marathon': '2019-10-13', 'New York City Marathon': '2019-11-03', 'Tokyo Marathon': '2019-03-03'},
    2021: {'Boston Marathon': '2021-10-11', 'London Marathon': '2021-10-03', 'Berlin Marathon': '2021-09-26', 'Chicago Marathon': '2021-10-10', 'New York City Marathon': '2021-11-07', 'Tokyo Marathon': '2021-03-06'},
    2022: {'Boston Marathon': '2022-04-18', 'London Marathon': '2022-10-02', 'Berlin Marathon': '2022-09-25', 'Chicago Marathon': '2022-10-09', 'New York City Marathon': '2022-11-06', 'Tokyo Marathon': '2022-03-06'},
    2023: {'Boston Marathon': '2023-04-17', 'London Marathon': '2023-04-23', 'Berlin Marathon': '2023-09-24', 'Chicago Marathon': '2023-10-08', 'New York City Marathon': '2023-11-05', 'Tokyo Marathon': '2023-03-05'},
    2024: {'Boston Marathon': '2024-04-15', 'London Marathon': '2024-04-21', 'Berlin Marathon': '2024-09-29', 'Chicago Marathon': '2024-10-13', 'New York City Marathon': '2024-11-03', 'Tokyo Marathon': '2024-03-03'},
    2025: {'Boston Marathon': '2025-04-21', 'London Marathon': '2025-04-27', 'Tokyo Marathon': '2025-03-02'}
}


def load_marathon_data(csv_path='marathon_results_top10.csv'):
    """Load marathon data and extract unique race dates"""
    df = pd.read_csv(csv_path)
    unique_races = df[['Marathon', 'Year']].drop_duplicates()
    
    race_info = []
    for _, row in unique_races.iterrows():
        marathon = row['Marathon']
        year = row['Year']
        
        if marathon in MARATHON_LOCATIONS and year in MARATHON_DATES:
            race_date = MARATHON_DATES[year].get(marathon)
            if race_date:
                info = {
                    'marathon': marathon,
                    'year': year,
                    'race_date': race_date,
                    'city': MARATHON_LOCATIONS[marathon]['city'],
                    'country': MARATHON_LOCATIONS[marathon]['country'],
                    'latitude': MARATHON_LOCATIONS[marathon]['lat'],
                    'longitude': MARATHON_LOCATIONS[marathon]['lon']
                }
                race_info.append(info)
    
    return pd.DataFrame(race_info)


def get_weather_data(lat, lon, date):
    """
    Get weather data from Open-Meteo API for a specific date
    Returns hourly data for marathon hours (6 AM - 4 PM)
    """
    # Open-Meteo Historical Weather API
    url = "https://archive-api.open-meteo.com/v1/archive"
    
    params = {
        'latitude': lat,
        'longitude': lon,
        'start_date': date,
        'end_date': date,
        'hourly': [
            'temperature_2m',
            'relative_humidity_2m',
            'precipitation',
            'wind_speed_10m',
            'wind_direction_10m',
            'surface_pressure'
        ],
        'timezone': 'auto'
    }
    
    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
        data = response.json()
        
        if 'hourly' not in data:
            return None
        
        hourly = data['hourly']
        
        # Extract marathon hours (typically 6 AM - 4 PM, indices 6-16)
        # Average the conditions during typical marathon running hours
        marathon_hours_start = 6
        marathon_hours_end = 16
        
        def avg_values(values, start_idx=marathon_hours_start, end_idx=marathon_hours_end):
            """Average values during marathon hours, handling None values"""
            valid_values = [v for v in values[start_idx:end_idx+1] if v is not None]
            return sum(valid_values) / len(valid_values) if valid_values else None
        
        result = {
            'temperature_c': avg_values(hourly['temperature_2m']),
            'humidity_percent': avg_values(hourly['relative_humidity_2m']),
            'wind_speed_ms': avg_values(hourly['wind_speed_10m']),
            'wind_direction_deg': avg_values(hourly['wind_direction_10m']),
            'precipitation_mm': sum([p for p in hourly['precipitation'][marathon_hours_start:marathon_hours_end+1] if p is not None]),
            'pressure_hpa': avg_values(hourly['surface_pressure'])
        }
        
        return result
        
    except Exception as e:
        print(f"  Error fetching weather: {e}")
        return None


def main():
    """Main execution"""
    
    print("=" * 80)
    print("World Marathon Majors - Weather Data Collection (Open-Meteo)")
    print("=" * 80)
    print()
    
    # Load marathon data
    print("Loading marathon race data...")
    races_df = load_marathon_data()
    print(f"Found {len(races_df)} marathon races with dates")
    print()
    
    print("=" * 80)
    print("Starting weather data collection...")
    print("Source: Open-Meteo Historical Weather API (free, no key needed)")
    print("Variables: temperature, humidity, wind speed/direction, precipitation, pressure")
    print("=" * 80)
    print()
    
    all_weather = []
    
    for idx, race in races_df.iterrows():
        marathon = race['marathon']
        year = race['year']
        race_date = race['race_date']
        lat = race['latitude']
        lon = race['longitude']
        city = race['city']
        
        print(f"[{idx+1}/{len(races_df)}] {marathon} {year}")
        print(f"  Location: {city} ({lat:.4f}, {lon:.4f})")
        print(f"  Race date: {race_date}")
        
        # Get weather data
        weather = get_weather_data(lat, lon, race_date)
        
        if weather:
            result = {
                'marathon': marathon,
                'year': year,
                'city': city,
                'country': race['country'],
                'race_date': race_date,
                **weather
            }
            all_weather.append(result)
            print(f"  ✓ Temperature: {weather['temperature_c']:.1f}°C, "
                  f"Humidity: {weather['humidity_percent']:.0f}%, "
                  f"Wind: {weather['wind_speed_ms']:.1f} m/s")
        else:
            print(f"  ✗ Failed to collect weather data")
        
        # Be polite to API
        time.sleep(0.5)
        print()
    
    # Save results
    if all_weather:
        weather_df = pd.DataFrame(all_weather)
        output_file = 'marathon_weather_data.csv'
        weather_df.to_csv(output_file, index=False)
        
        print("=" * 80)
        print("COLLECTION COMPLETE")
        print("=" * 80)
        print(f"✓ Collected weather data for {len(all_weather)} races")
        print(f"✓ Saved to: {output_file}")
        print()
        print("Sample data:")
        print(weather_df.head(10).to_string())
        print()
    else:
        print("✗ No weather data collected")
        sys.exit(1)


if __name__ == '__main__':
    main()
