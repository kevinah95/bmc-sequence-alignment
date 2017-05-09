import atexit
from time import time
from datetime import timedelta

def secondsToStr(t):
    return str(timedelta(seconds=t))

line = "="*40
def log(s, elapsed=None):
    print(line)
    print(secondsToStr(time()), '-', s)
    if elapsed:
        print("Elapsed time:", elapsed)
    print(line)
    print()

def endlog():
    end = time()
    elapsed = end-start
    log("End Program", secondsToStr(elapsed))

def now():
    return secondsToStr(time())

start = time()
atexit.register(endlog)
log("Start Program")



endlog()#Se llama a la función
