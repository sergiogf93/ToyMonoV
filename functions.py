# Here the general functions will be defined

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    ERROR = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def printError(message):
    print(bcolors.ERROR + message + bcolors.ENDC)


def Get2DBin(hist, x, y):
    xaxis = hist.GetXaxis()
    yaxis = hist.GetYaxis()
    binx = xaxis.FindBin(x)
    biny = yaxis.FindBin(y)
    return hist.GetBin(binx,biny)
