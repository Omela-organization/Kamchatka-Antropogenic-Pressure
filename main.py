import schedule

from src.anropogenic_pressure import run

schedule.every(1).day.do(run)
while True:
    schedule.run_pending()
