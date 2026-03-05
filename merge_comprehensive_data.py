"""
Merge comprehensive marathon results with pollution and weather data
"""
import pandas as pd
import numpy as np

print("="*80)
print("MERGING COMPREHENSIVE DATA WITH POLLUTION AND WEATHER")
print("="*80)
print()

# Load data
print("Loading data...")
runners = pd.read_csv('marathon_comprehensive_api.csv', low_memory=False)
pollution = pd.read_csv('marathon_pollution_data.csv')
weather = pd.read_csv('marathon_weather_data.csv')

print(f"✓ Runners: {len(runners):,} rows")
print(f"✓ Pollution: {len(pollution):,} rows")
print(f"✓ Weather: {len(weather):,} rows")
print()

# Check pollution data structure
print("Pollution data columns:")
print(pollution.columns.tolist())
print()

# Prepare merge keys
print("Preparing merge keys...")

# Standardize marathon names
name_mapping = {
    'Boston Marathon': 'Boston',
    'London Marathon': 'London',
    'Chicago Marathon': 'Chicago',
    'New York City Marathon': 'NYC',
    'Berlin Marathon': 'Berlin',
    'Tokyo Marathon': 'Tokyo'
}

runners['Marathon'] = runners['marathon'].map(name_mapping)
pollution['Marathon'] = pollution['marathon'].map(name_mapping)
weather['Marathon'] = weather['marathon'].map(name_mapping)

# Filter pollution to only race day (not pre-race days)
if 'date_type' in pollution.columns:
    pollution = pollution[pollution['date_type'] == 'race_day'].copy()
    print(f"Filtered pollution to race_day only: {len(pollution)} rows")

print(f"Unique marathons in runners data: {sorted(runners['Marathon'].dropna().unique())}")
print(f"Unique marathons in pollution data: {sorted(pollution['Marathon'].dropna().unique())}")
print()

# Rename year column if needed
if 'year' in runners.columns:
    runners['Year'] = runners['year']
if 'year' in pollution.columns:
    pollution['Year'] = pollution['year']
if 'year' in weather.columns:
    weather['Year'] = weather['year']

# Merge with pollution
print("Merging with pollution data...")
merged = runners.merge(pollution[['Marathon', 'Year', 'NO2', 'PM2.5', 'PM10']],
                       on=['Marathon', 'Year'], how='left', suffixes=('', '_race_day'))
print(f"✓ After pollution merge: {len(merged):,} rows")

# Rename pollution columns to match expected names
if 'PM2.5' in merged.columns:
    merged.rename(columns={'PM2.5': 'PM2.5_race_day', 'PM10': 'PM10_race_day', 'NO2': 'NO2_race_day'}, inplace=True)

# Check match rate
matched_pollution = merged['PM2.5_race_day'].notna().sum()
print(f"  Matched pollution: {matched_pollution:,} ({matched_pollution/len(merged)*100:.1f}%)")
print()

# Merge with weather
print("Merging with weather data...")
weather_cols = [c for c in weather.columns if c not in ['Marathon', 'Year', 'marathon', 'year']]
merged = merged.merge(weather[['Marathon', 'Year'] + weather_cols], on=['Marathon', 'Year'], how='left')
print(f"✓ After weather merge: {len(merged):,} rows")

# Standardize weather column names
weather_rename = {}
for old, new in [('temp_2m_mean', 'temp_race_day'),
                 ('temperature_c', 'temp_race_day'),
                 ('relative_humidity_2m_mean', 'humidity_race_day'),
                 ('humidity_percent', 'humidity_race_day'),
                 ('wind_speed_10m_max', 'wind_speed_race_day'),
                 ('wind_speed_ms', 'wind_speed_race_day')]:
    if old in merged.columns and new not in merged.columns:
        weather_rename[old] = new

if weather_rename:
    merged.rename(columns=weather_rename, inplace=True)
    print(f"  Renamed weather columns: {weather_rename}")

if 'temp_race_day' in merged.columns:
    matched_weather = merged['temp_race_day'].notna().sum()
    print(f"  Matched weather: {matched_weather:,} ({matched_weather/len(merged)*100:.1f}%)")
else:
    print("  Warning: No temperature column found")
print()

# Convert time to minutes
print("Converting finish times to minutes...")

def parse_time(t):
    """Convert HH:MM:SS to minutes"""
    if pd.isna(t) or t is None:
        return np.nan
    try:
        parts = str(t).split(':')
        if len(parts) == 3:
            h, m, s = parts
            return int(h) * 60 + int(m) + float(s) / 60
        elif len(parts) == 2:
            m, s = parts
            return int(m) + float(s) / 60
        else:
            return float(t)
    except:
        return np.nan

merged['time_minutes'] = merged['final_time'].apply(parse_time)

valid_times = merged['time_minutes'].notna().sum()
print(f"✓ Valid finish times: {valid_times:,} ({valid_times/len(merged)*100:.1f}%)")
print()

# Summary statistics
print("="*80)
print("FINAL DATASET SUMMARY")
print("="*80)
print(f"Total rows: {len(merged):,}")
print()

print("By Marathon:")
print(merged.groupby('Marathon').size().sort_values(ascending=False))
print()

print("By Year:")
print(merged.groupby('Year').size().sort_index())
print()

print("Data completeness:")
print(f"  Finish time: {merged['time_minutes'].notna().sum():,} ({merged['time_minutes'].notna().sum()/len(merged)*100:.1f}%)")
print(f"  Age: {merged['age'].notna().sum():,} ({merged['age'].notna().sum()/len(merged)*100:.1f}%)")
print(f"  Sex: {merged['sex'].notna().sum():,} ({merged['sex'].notna().sum()/len(merged)*100:.1f}%)")
print(f"  PM2.5: {merged['PM2.5_race_day'].notna().sum():,} ({merged['PM2.5_race_day'].notna().sum()/len(merged)*100:.1f}%)")
if 'temp_race_day' in merged.columns:
    print(f"  Temperature: {merged['temp_race_day'].notna().sum():,} ({merged['temp_race_day'].notna().sum()/len(merged)*100:.1f}%)")
else:
    print(f"  Temperature: No data (column not found)")
print()

# Save
output_file = 'marathon_comprehensive_with_pollution.csv'
print(f"Saving to {output_file}...")
merged.to_csv(output_file, index=False)

print(f"✓ Saved: {output_file}")
print(f"  Size: {len(merged):,} rows")

# Also save a summary
agg_dict = {
    'time_minutes': ['count', 'mean', 'median', 'std'],
    'PM2.5_race_day': 'first',
    'PM10_race_day': 'first',
    'NO2_race_day': 'first'
}
if 'temp_race_day' in merged.columns:
    agg_dict['temp_race_day'] = 'first'
if 'humidity_race_day' in merged.columns:
    agg_dict['humidity_race_day'] = 'first'

summary = merged.groupby(['Marathon', 'Year']).agg(agg_dict).reset_index()

summary.to_csv('marathon_race_summary.csv', index=False)
print(f"✓ Saved: marathon_race_summary.csv (race-level summary)")

print()
print("="*80)
print("MERGE COMPLETE")
print("="*80)
