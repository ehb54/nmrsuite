import time
class Stopwatch:
    def __init__(self):
        self.paused = False
        self.starttime = float('nan')
    def start (self):
        self.starttime = time.time()
    def pause (self): # time stops running when paused
        if (self.paused == False):
            self.pausestart = time.time()
    def resume (self): 
        if (self.paused == True):
            self.paused = False
            self.starttime += (time.time() - self.pausestart)
    def get(self):
        print ("Elapsed Time:", (time.time() - self.starttime), "seconds")

'''
Sample:
s = Stopwatch()
s.start()
... (2 seconds of code) ...
s.pause()
... (1 second of code) ...
s.resume()
... (3 seconds of code) ...
s.get() - prints "Elapsed Time: 5 seconds"
'''