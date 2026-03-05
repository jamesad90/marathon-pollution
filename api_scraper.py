"""
Fast API-based Marathon Scraper
Uses marathonguide.com's backend API (runzy.com) directly
MUCH faster than Selenium!
"""

import requests
import pandas as pd
import time
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class FastAPIMarathonScraper:
    """Fast scraper using runzy.com API"""

    def __init__(self):
        self.api_base = "https://back.runzy.com/mg/event-results"
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
        })

        # Marathon base slugs and their year offset patterns
        # Each marathon has a linear pattern from 2010-2024
        self.marathons = {
            'Boston Marathon': ('boston-marathon', 1999),
            'London Marathon': ('london-marathon', 1999),  # Will verify base
            'Chicago Marathon': ('chicago-marathon', 1999),  # Will verify base
            'Berlin Marathon': ('berlin-marathon', 1999),  # Will verify base
            'Tokyo Marathon': ('tokyo-marathon', 1999),  # Will verify base
            'New York City Marathon': ('new-york-city-marathon', 1998),
        }

        # Special case slugs for 2025 (not following linear pattern)
        self.slug_2025 = {
            'Boston Marathon': 'boston-marathon-22',
            'London Marathon': 'london-marathon-21',
            'Chicago Marathon': 'chicago-marathon-21',  # Need to verify
            'Berlin Marathon': 'berlin-marathon-21',  # Need to verify
            'Tokyo Marathon': 'tokyo-marathon-21',  # Need to verify
            'New York City Marathon': 'new-york-city-marathon-21',  # Need to verify
        }

    def get_event_slug(self, marathon_name: str, year: int) -> str:
        """
        Generate the correct event slug with year-based number suffix.

        Pattern for 2010-2024: {base-slug}-{year - base_year}
        Exception: 2023 has no number suffix
        Exception: 2025 uses special mapping (not linear)

        Examples:
        - Boston 2024 = boston-marathon-25 (2024-1999=25)
        - Boston 2023 = boston-marathon (no suffix)
        - Boston 2025 = boston-marathon-22 (special case)
        """
        if marathon_name not in self.marathons:
            return None

        # Special case for 2025
        if year == 2025:
            return self.slug_2025.get(marathon_name)

        # Special case for 2023 - no number suffix
        if year == 2023:
            base_slug, _ = self.marathons[marathon_name]
            return base_slug

        # Normal case for 2010-2024 (except 2023)
        base_slug, base_year = self.marathons[marathon_name]
        number = year - base_year
        return f"{base_slug}-{number}"

    def fetch_page(self, event_slug: str, year: int, page: int, limit: int = 100, max_retries: int = 3):
        """Fetch a single page from API with retry logic for rate limiting"""

        url = f"{self.api_base}/{event_slug}/"

        params = {
            'subevent': 'all',
            'gender': 'all',
            'age_group': 'all',
            'page': page,
            'limit': limit,
            'order_by': 'over_all_place',
            'order_dir': 'asc',
            'year': year
        }

        for attempt in range(max_retries):
            try:
                response = self.session.get(url, params=params, timeout=15)

                if response.status_code == 200:
                    data = response.json()
                    return {
                        'success': True,
                        'results': data.get('results', []),
                        'pagination': data.get('pagination', {}),
                        'summary': data.get('summary', {})
                    }
                elif response.status_code == 429:
                    # Rate limited - wait and retry
                    wait_time = (attempt + 1) * 2  # Exponential backoff: 2, 4, 6 seconds
                    logger.warning(f"Rate limited on page {page}, waiting {wait_time}s...")
                    time.sleep(wait_time)
                    continue
                else:
                    logger.warning(f"Page {page}: HTTP {response.status_code}")
                    return {'success': False, 'results': [], 'pagination': {}, 'summary': {}}

            except Exception as e:
                if attempt < max_retries - 1:
                    time.sleep(1)
                    continue
                logger.warning(f"Error fetching page {page}: {e}")
                return {'success': False, 'results': [], 'pagination': {}, 'summary': {}}

        # All retries exhausted
        return {'success': False, 'results': [], 'pagination': {}, 'summary': {}}

    def scrape_marathon(self, marathon_name: str, year: int, max_workers: int = 1):
        """
        Scrape entire marathon using API calls

        marathon_name: e.g. 'Boston Marathon'
        year: e.g. 2024
        max_workers: parallel requests (default 1 to avoid rate limiting)
        """

        event_slug = self.get_event_slug(marathon_name, year)
        if not event_slug:
            logger.error(f"Unknown marathon: {marathon_name}")
            return pd.DataFrame()

        logger.info(f"Scraping {marathon_name} {year}")

        # Get first page to determine total
        first_page = self.fetch_page(event_slug, year, page=1, limit=100)

        if not first_page['success']:
            logger.warning(f"{marathon_name} {year}: No data")
            return pd.DataFrame()

        pagination = first_page['pagination']
        total = pagination.get('total', 0)
        last_page = pagination.get('last_page', 0)

        logger.info(f"  Total finishers: {total:,}, Pages: {last_page}")

        if total == 0:
            return pd.DataFrame()

        # Collect all results
        all_results = []
        all_results.extend(first_page['results'])

        # Add delay after first page if sequential
        if max_workers == 1:
            time.sleep(0.5)

        # Fetch remaining pages
        if last_page > 1:
            pages_to_fetch = list(range(2, last_page + 1))

            if max_workers == 1:
                # Sequential fetch with delay to avoid rate limiting
                for page in pages_to_fetch:
                    try:
                        page_data = self.fetch_page(event_slug, year, page, 100)
                        if page_data['success']:
                            all_results.extend(page_data['results'])

                        if page % 50 == 0:
                            logger.info(f"    Progress: {page}/{last_page} pages")

                        # Delay to avoid rate limiting
                        time.sleep(0.5)

                    except Exception as e:
                        logger.error(f"Error processing page {page}: {e}")
            else:
                # Parallel fetch (may hit rate limits)
                with ThreadPoolExecutor(max_workers=max_workers) as executor:
                    futures = {
                        executor.submit(self.fetch_page, event_slug, year, page, 100): page
                        for page in pages_to_fetch
                    }

                    for future in as_completed(futures):
                        page_num = futures[future]
                        try:
                            page_data = future.result()
                            if page_data['success']:
                                all_results.extend(page_data['results'])

                            if page_num % 50 == 0:
                                logger.info(f"    Progress: {page_num}/{last_page} pages")

                        except Exception as e:
                            logger.error(f"Error processing page {page_num}: {e}")

        # Convert to DataFrame
        df = pd.DataFrame(all_results)

        if not df.empty:
            df['year'] = year
            df['marathon'] = marathon_name

        logger.info(f"  ✓ Collected {len(df):,} results")

        return df

    def scrape_all_marathons(self, years: list, output_file: str = 'marathon_api_results.csv',
                           max_workers: int = 3):
        """
        Scrape all marathons for given years using threading

        years: list like [2020, 2021, 2022, 2023, 2024]
        max_workers: number of marathons to scrape in parallel
        """

        all_data = []
        data_lock = Lock()

        # Create tasks
        tasks = []
        for marathon in self.marathons.keys():
            for year in years:
                tasks.append((marathon, year))

        print(f"="*80)
        print(f"FAST API SCRAPER")
        print(f"="*80)
        print(f"Total tasks: {len(tasks)}")
        print(f"Parallel scrapers: {max_workers}")
        print()

        def scrape_task(marathon, year):
            """Thread worker"""
            try:
                df = self.scrape_marathon(marathon, year)

                if not df.empty:
                    with data_lock:
                        all_data.append(df)

                        # Save checkpoint
                        if all_data:
                            combined = pd.concat(all_data, ignore_index=True)
                            combined.to_csv('checkpoint_' + output_file, index=False)
                            logger.info(f"[Checkpoint] {len(combined):,} total results")

                return (marathon, year, len(df) if not df.empty else 0)

            except Exception as e:
                logger.error(f"Error in {marathon} {year}: {e}")
                return (marathon, year, 0)

        # Run in parallel
        completed = 0
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(scrape_task, m, y): (m, y) for m, y in tasks}

            for future in as_completed(futures):
                marathon, year, count = future.result()
                completed += 1
                print(f"[Progress] {completed}/{len(tasks)}: {marathon} {year} ({count:,} results)")

        # Final save
        if all_data:
            combined = pd.concat(all_data, ignore_index=True)
            combined.to_csv(output_file, index=False)

            print()
            print(f"="*80)
            print(f"COMPLETE!")
            print(f"="*80)
            print(f"Total results: {len(combined):,}")
            print(f"Saved to: {output_file}")
            print()
            print("By marathon:")
            print(combined.groupby('marathon').size().sort_values(ascending=False))
            print()
            print("By year:")
            print(combined.groupby('year').size().sort_index())

            return combined

        return pd.DataFrame()


def test():
    """Quick test"""
    print("="*80)
    print("TESTING API SCRAPER - Boston Marathon 2024")
    print("="*80)
    print()

    scraper = FastAPIMarathonScraper()
    df = scraper.scrape_marathon('Boston Marathon', 2024)

    if not df.empty:
        print(f"\n✓ Success! {len(df):,} results")
        print()
        print("Sample:")
        print(df.head(20))
        print()
        print(f"Columns: {df.columns.tolist()}")

        df.to_csv('test_api_scraper.csv', index=False)
        print("\n✓ Saved to test_api_scraper.csv")
    else:
        print("\n✗ No results - may need to find correct event slug")


def main():
    """Full scrape"""
    scraper = FastAPIMarathonScraper()

    # Scrape 2010-2024 (skip 2025 for now)
    years = list(range(2010, 2025))  # 2010-2024

    # Use fewer workers to avoid rate limiting
    scraper.scrape_all_marathons(years, output_file='marathon_comprehensive_api.csv', max_workers=2)


if __name__ == '__main__':
    # Run full scrape
    main()
