import time

# timeLabels will be set in the following manner: timeLabels = { "timeLabel" : {"elapsedTime" : 0.0, "initialTime" : 0.0} }
class Timer:
	def __init__(self, timeStep):
		self.timeStep = timeStep
		self.timeLabels = dict()
		self.currentTime = 0.0

	def start(self, timeLabel):
		self.timeLabels[timeLabel] = { "elapsedTime" : 0.0 , "initialTime" : time.time() }

	def stop(self, timeLabel):
		self.timeLabels[timeLabel]["elapsedTime"] += time.time() - self.timeLabels[timeLabel]["initialTime"]
		self.timeLabels[timeLabel]["initialTime"] = 0.0

	def incrementTime(self):
		self.currentTime += self.timeStep

	def getCurrentTime(self):
		return self.currentTime
