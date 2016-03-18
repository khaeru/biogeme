#!/usr/bin/python3

import os, time, threading
import subprocess
import inspect
import html.parser
import zipfile
import datetime
import configparser
import sys
import webbrowser
import re
import socket

onWindows = False
if sys.platform == "nt" or sys.platform == "win32":
    onWindows = True



#trying to import cherrypy if not available get it via pip
#if pip not available, get pip first
try:
    import cherrypy
except ImportError:
    getpip = ["python3","get-pip.py"]
    getcpy = ["python3","-m","pip","install","cherrypy"]
    if onWindows:
       getpip[0] = "python"
       getcpy[0] = "python"
    else:
        getpip = ["sudo"] + getpip
        getcpy = ["sudo"] + getcpy

    print("CherryPy is missing, trying to install")
    #stderr is not redirected to subprocess.PIPE so the user can see the errors
    print(getpip)
    subprocess.call(getpip, stdout=subprocess.PIPE, shell=onWindows)
    print("pip installed")
    print(getcpy)
    subprocess.call(getcpy,stdout=subprocess.PIPE, shell=onWindows)
    print("cherrypy installed")
    import cherrypy




#resources location
localDir = os.path.dirname(__file__)
absDir = os.path.join(os.getcwd(), localDir)
BS_FOLDER =  os.path.join(os.getcwd(), "bootstrap")
JBS_FOLDER =  os.path.join(os.getcwd(), "jasny-bootstrap")
CSS_FOLDER =  os.path.join(os.getcwd(), "css")
DOXY_FOLDER = os.path.join(os.getcwd(), "html/doxygen/html")
IMG_FOLDER =  os.path.join(os.getcwd(), "img")
ICONE = os.path.join(os.getcwd(), "img/B.ico")



#reading config from file
config = configparser.SafeConfigParser()
config.read("biogemeWebUI.conf")
configServer = configparser.SafeConfigParser()
configServer.read("server.conf")

ef = configServer["global"]["log.error_file"]
if not os.path.exists(ef):
    # 'touch' the file if it doesnt exists
    print("create log file " + ef[1:-1])
    open(os.path.join(os.getcwd(), ef[1:-1]),"w").close()


DEFAULT_SORT = config["biogemeWebUI"]["defaultSort"]
print("\n****** Using '%s' to sort ******" % DEFAULT_SORT)

df = config["biogemeWebUI"]["dataFolder"]
DATA_FOLDER = os.path.join(os.getcwd(), df)
print("****** Using '%s' as data folder ******" % DATA_FOLDER)
if not os.path.isdir(DATA_FOLDER):
    print("- Making directory")
    os.makedirs(DATA_FOLDER)


#the webui should not be SHUTDOWNABLE if running on a distant server
SHUTDOWNABLE = False
if config["biogemeWebUI"]["shutdown_allowed"] == "True":
    SHUTDOWNABLE = True
print("****** Shutdown enabled: %s ******\n" % SHUTDOWNABLE)



#web app that will be loaded by cherrypy
class BiogemeWebUI(object):
    #process that will run the commands
    p = None
    #if the webui is busy
    busy = False
    lastModelUsed = ""
    lastDataUsed = ""
    
    


    def upload(self, myFile):
        #TODO
        #check if file already exists

        #check if data folder is present
        if not os.path.exists(DATA_FOLDER):
            os.makedirs(DATA_FOLDER)

        #create file
        if myFile.filename != "":

            newFile = open(os.path.join(DATA_FOLDER, myFile.filename), "wb")

            #read data by chunks
            size = 0
            while True:
                data = myFile.file.read(8192)
                if not data:
                    break
                size += len(data)
                #writing to file
                newFile.write(data)

            return self.prepareAlert( myFile.filename + " has been uploaded" )
        else:
            return self.prepareAlert( "Please select a file first", "danger")

    def test(self):
        cherrypy.engine.restart()
        return "Hello, you tried to restart the cherrypy engine."
    test.exposed = True


    #get all files in the DATA_FOLDER with name and date
    #can be filter by extension and/or basename
    def listDataFolder(self, extension, basename=""):
        fileList = []
        global DATA_FOLDER
        print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" + DATA_FOLDER)
        for f in os.listdir(DATA_FOLDER):
            if f.endswith(extension) and f.startswith(basename):
                filePath = os.path.join(DATA_FOLDER, f)
                if not os.path.isdir(filePath):
                    fileList.append([f, time.strftime('%a %d.%m.%y %H:%M:%S', time.gmtime(os.path.getmtime(filePath)))])
        return fileList     


    def callBiogeme(self, cmd, modelFile, dataFile):
        command = "echo"
        global DATA_FOLDER
        
        if not self.busy:
            commands = []
            #ping is only for debug purposes
            if cmd.startswith("ping"):
                count = cmd[4:]
                try:
                    int(count)
                except ValueError:
                    #some tried something weird
                    return

                commands = ["ping", "-c", count, "127.0.0.1"]
                if onWindows:
                    commands [1] = "-n"
                self.lastModelUsed="ping"
                self.lastDataUsed="ping"

            else:
                if cmd == "0":
                    command = "echo"
                elif cmd == "1":
                    command = "pythonbiogeme"
                elif cmd == "2":
                    command = "biogeme"
                elif cmd == "3":
                    command = "biosim"

                print("=========================== Running ========================")
                print(command + " " + modelFile + " " + dataFile)

                commands = [command, modelFile, dataFile]
                self.lastModelUsed=modelFile
                self.lastDataUsed=dataFile

            self.p = subprocess.Popen(commands, cwd=DATA_FOLDER, stdout=subprocess.PIPE,shell=onWindows)

            self.busy = True

            return True
        else:
            #already running
            print("already running")
            return False


    #used to kill a running command
    #could be more secure
    def stopBiogeme(self, pw=""):
        if pw == "password":
            self.p.kill()
            self.busy = False
            self.lastModelUsed = ""
            self.lastDataUsed = ""
    stopBiogeme.exposed = True
    

    #read a bit of the console output of the command
    def readBuffer(self, nbLineToRead=0, ms=1000):
        cpt = 0
        rep = "event: console-output\ndata: "
        start = int(round(time.time()*1000))
        while start + ms >int(round(time.time()*1000)):
            line = self.p.stdout.readline()
            
            print(line)
            #line != b'' and
            if line != b'' and (cpt < nbLineToRead or nbLineToRead == 0):
                if os.name == "nt":
                    decoded = line.decode("windows-1252")
                    if decoded.find("\r\n") != -1:
                        rep += decoded.replace("\r\n", "<br />")
                    else:
                        rep += decoded.replace("\n", "<br />")
                elif os.name == "posix":
                    rep += line.decode("utf-8").replace("\n", "<br />")
                cpt += 1
            else:
                print("**************************************************** break")
                break
        return rep


    #prepare a link to the results
    def linkToResults(self):
        if self.lastModelUsed != "":
            return "<a href='results?sort=dateDes&filter="+self.lastModelUsed+"'>See results for "+self.lastModelUsed+"</a>"
        else:
            return "<a href='results'>See results</a>"

    
    #prepare html tables for the run pages
    def getDragableTable(self, extension,basename="",sort=DEFAULT_SORT,baseLink="",checkboxes=False):
        ext_len = len (extension) 
        
        table = """
    <table class="ulsd">
        <tr>"""
        if checkboxes:
            table+="""
            <th><i id="checkAll" class="glyphicon glyphicon-ok-circle" value="true"></i></th>"""


        table += """
            <th>Name <a href="""+ baseLink +"""?sort=nameAsc><span class="glyphicon glyphicon-chevron-up"></span></a><a href="""+ baseLink +"""?sort=nameDes><span class="glyphicon glyphicon-chevron-down"></a></th>
            <th>Last modified <a href="""+ baseLink +"""?sort=dateAsc><span class="glyphicon glyphicon-chevron-up"></a><a href="""+ baseLink +"""?sort=dateDes><span class="glyphicon glyphicon-chevron-down"></a></th>
        </tr>
    </table>"""



        table += """
<div class="listHolder">
    <table class="ulsd"  id="filesTable">
        """
        files = self.listDataFolder(extension,basename)

        if sort == "nameAsc":
            files.sort(key=lambda f: f[0])
        elif sort == "dateAsc":
            files.sort(key=lambda f: f[1])
        elif sort == "nameDes":
            files.sort(key=lambda f: f[0],reverse=True)
        elif sort == "dateDes":
            files.sort(key=lambda f: f[1],reverse=True)

        i = 1

        for element in files:
            table += """
        <tr>"""
            if checkboxes:
                table+="""<td><input type="checkbox" value='"""+ element[0] +  """' id='cb"""+ str(i) + """'></td>"""
                i+=1

            table+="""<td><a href="data/"""
            table += element[0]
            table += """" draggable="true" ondragstart="drag(event)">"""
            if(ext_len > 0):
                table += element[0][:-ext_len] + """</a></td>"""
            else:
                table += element[0] + """</a></td>"""
            table += """<td>""" + element[1] + """</td></tr>"""
        table += """
    </table>
</div>"""
        return table

    def egg(self):
        f = open("html/help.html")
        page = f.read()

        d = page.index("<deb></deb>")
        e = page.index("<end></end>") + len("<end></end>")
        page = page[:d] + "<marquee>°</marquee>" + page[e:]
     
        footer = self.getFooter()
        page = page.replace("DYNAMIC_FOOTER", footer)
        navBar = self.getNavBar()
        page = page.replace("DYNAMIC_NAV_BAR",navBar)
        return page
    egg.exposed = True

    #prepare the footer
    def getFooter(self):
        out = """
        <!--    Footer      -->
        <div class='altfoot'>
            <div class="col-md-4 alignLeft" id="FOOTER_LEFT"></div>
            <div class="col-md-4 alignCenter" id="FOOTER_CENTER"></div>
            <div class="col-md-4 alignRight" id="FOOTER_RIGHT"></div>
        </div> """

        return out

    #prepare the navBar, possibility to highlight item and subitem
    def getNavBar(self,itemNumber=0, subitemNumber=0):
        navBar = """
        <div class="navbar navbar-default navbar-fixed-top" role="navigation">
            <div class="container">
                <div class="navbar-header">
                    <button type="button" class="navbar-toggle" data-toggle="collapse" data-target=".navbar-collapse">
                    <span class="sr-only">Toggle navigation</span>
                    <span class="icon-bar"></span>
                    <span class="icon-bar"></span>
                    <span class="icon-bar"></span>
                    </button>
                    <a class="navbar-brand" href="home">Biogeme web ui</a>
                </div>
                <div class="collapse navbar-collapse">
                    <ul class="nav navbar-nav">
                        <li><a href="files">Files</a></li>
                        <li class="dropdown"><a href="#" class="dropdown-toggle" data-toggle="dropdown">Run <b class="caret"></b></a>
                            <ul class="dropdown-menu">
                            <li><a href="runpb">Pythonbiogeme</a></li>
                            <li><a href="runb">Biogeme</a></li>
                            <li><a href="runbs">Biosim</a></li>
                            </ul>
                        </li>
                        <li><a href="console">Console</a></li>
                        <li><a href="results">Results</a></li>
                    </ul>
                    <ul class="nav navbar-nav navbar-right">
                        <li><a href="doc"><span class="glyphicon glyphicon-book" title="Documentation"></a></li>
                        <li><a href="help"><span class="glyphicon glyphicon-question-sign" title="Help"></a></li>"""
        if SHUTDOWNABLE:
            navBar += """
                        <li><a href="settings"><span class="glyphicon glyphicon-cog" title="Settings"></span></a></li>
                        <li><a href="shutdown"><span class="glyphicon glyphicon-off" title="Shut the webui down"></a></li>"""
                    
        navBar+="""
                    </ul>
                </div><!--/.nav-collapse -->
              </div>
            </div>
        """

        #ugly inner-methods to highlight (sub)items
        def nS(n, i):
            try:
                r = n.index("<li", i)
            except ValueError:
                r = 10000
            return r

        def nE(n, i):
            return n.index("</li", i)

        def findPos(n, item,start=0):
            level = 1
            pos = nS(n, start)
            i = {0:1, 1:1, 2:1, 3:1}

            while i[0] < item:
                s = nS(n, pos+1)
                e = nE(n, pos+1)

                if s < e:
                    level += 1
                    pos = s
                else:
                    i[level] = 1
                    level -= 1
                    pos = e

                i[level] = i[level] + 1

            return nS(n, pos)

        def activateItem(navBar, p):
            if navBar[p:p+4] == "<li>":
                navBar = navBar[:p+3] + ' class="active"' + navBar[p+3:] 
            else:
                try:
                    if navBar.index('class="',p) < navBar.index(">",p):
                        c = navBar.index('class="',p) + len('class="')
                        navBar = navBar[:c] + 'active ' + navBar[c:]
                except ValueError:
                    print("index of navbar probably to high")

            return navBar

        def activate(nav, item, subitem=0):
            r= nav
            if item > 0:
                p = findPos(r, item)
                r = activateItem(r, p)

                if subitem > 0:
                    s = findPos(r,subitem,p+1)
                    r = activateItem(r, s)
            return r


        return activate(navBar, itemNumber, subitemNumber)

    #prepare the zip file to be downloaded
    def prepareZip(self,files):
        global DATA_FOLDER
        now = datetime.datetime.now()
        name = "Files_"
        name +=str(now.date())
        name +="_"
        name +=str(now.hour) +"-" +str("%02i" % now.minute)
        name +=".zip"

        zp = os.path.join(DATA_FOLDER, "zip")
        if not os.path.isdir(zp):        
            os.makedirs(zp)
        zp = os.path.join(zp, name)
        zf = zipfile.ZipFile(zp,"w")
        for f in files:
            p = os.path.join(DATA_FOLDER, f)
            zf.write(p,f)
        zf.close()
        return "data/zip/" + name

    #tries to delete the files in the list
    def deleteFiles(self,files):
        global DATA_FOLDER
        success = True
        for f in files:
            p = os.path.join(DATA_FOLDER, f)
            try:
                os.remove(p)
            except FileNotFoundError:
                print(p + " was not found")
                success =False
        
        return success

    #prepare a bootstrap alert, acts like a notification after an action has succeeded or failed
    def prepareAlert(self, message, level="success"):
        return """
<div class="row marginTop">
    <div id="BS_ALERT" class="alert alert-"""+level+""" alert-dismissable">
        <button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>
        """+ message + """
    </div>
</div>"""


















    
#======================= SSE sources =========================================

    #sends console output and force to refresh at some rate
    def sendConsoleOutput(self):
        print("===================================================================")
        cherrypy.response.headers["Content-Type"] = "text/event-stream"
        cherrypy.response.headers["Cache-Control"] = "no-cache"

        refreshRate = 1000

        rep = "retry: "+ str(refreshRate) +"\n"
        if self.p == None:
            rep = "event: problem\ndata: "
            rep += "No process launched<br />\n\n"
        elif self.busy == False:
            #close connection
            rep = "event: end\ndata: "
            rep += self.linkToResults()
        elif self.p.poll() != None:
            #process done, reading the rest of output
            rep = self.readBuffer()
            self.busy = False
        else:
            #running, reading console output for a certin amount of time
            rep = self.readBuffer(ms=refreshRate)

        rep+="\n\n"
        return rep
    sendConsoleOutput.exposed = True


    #sends the status, free or busy. forces a refresh rate.
    def status(self):
        cherrypy.response.headers["Content-Type"] = "text/event-stream"
        cherrypy.response.headers["Cache-Control"] = "no-cache"
        
        refreshRate = 1000

        rep = "retry: "+ str(refreshRate) +"\n"
        rep += "event: status\n"
        if self.p == None or self.p.poll() != None:
            rep += "data: free\n\n"
        else:
            rep += "data: busy\n\n"
        return rep
    status.exposed = True















    
#======================= PAGES =========================================

    #home page, contains a frame containing biogeme.epfl.ch
    def home(self):
        f = open("html/home.html")
        page = f.read()
        
        footer = self.getFooter()
        page = page.replace("DYNAMIC_FOOTER", footer)
        navBar = self.getNavBar()
        page = page.replace("DYNAMIC_NAV_BAR", navBar)
        return page
    home.exposed = True

    #called when the root of the webui is called
    #links to the files page
    def index(self):
        return self.files()
    index.exposed = True


    #Displays all files in the working directory.
    #-sort can be used to define how to sort the files (nameAsc, nameDes, dateAsc, dateDes)
    #-filter allows to filter the files to display
    #-download can contain a list of files to zip and download
    #-delete can contain a list of files to remove from the working directory
    #-newFile can contain a file to add to the working directory
    def files(self,sort=DEFAULT_SORT, filter="",download="",delete="",newFile=None):
        f = open('html/files.html')
        page = f.read()


        msg = ""
        if newFile != None:
            msg += self.upload(newFile)
            


        if download != "":
            html_parser = html.parser.HTMLParser()
            unescaped = html_parser.unescape(download)
            files = unescaped.split(",")
            l = self.prepareZip(files)
            m = """<a href ='"""+l+"""' class="alert-link">Zip of files</a> should have been downloaded."""
            msg += self.prepareAlert(m)
            msg += """
            <iframe width="1" height="1" frameborder="0" src='"""+l+"""'></iframe>"""


        if delete != "":
            html_parser = html.parser.HTMLParser()
            unescaped = html_parser.unescape(delete)
            files = unescaped.split(",")
            delResult = self.deleteFiles(files)
            l = len(files)
            maxFiles = 5
            m=""
            if l == 1:
                m = files[0] + " has been deleted."
            elif l < maxFiles:
                for i in range(0,l):
                    m += files[i]
                    if i == l - 2:
                        m+= " and "
                    elif i < l -2:
                        m += ", "
                m+= " have been deleted"
            else:
                m = str(l) + " files have been deleted."

            if delResult:
                msg += self.prepareAlert(m)
            else:
                msg += self.prepareAlert(m,"warning")

            "" #files have been deleted"    + str(files)

        #filter the files to display
        ext = ""
        if filter.startswith("*"):
            ext = filter[1:]
            filter = ""
        allFiles = self.getDragableTable(ext,filter,sort, inspect.stack()[0][3],True)
        navBar = self.getNavBar(1)
        footer = self.getFooter()

        page = page.replace("ALL_FILES", allFiles)
        page = page.replace("MESSAGE", msg)
        page = page.replace("DYNAMIC_FOOTER", footer)
        page = page.replace("DYNAMIC_NAV_BAR", navBar)
        return page
    files.exposed = True

    #page to choose the model and data to launch with pythonbiogeme
    def runpb(self,sort=DEFAULT_SORT,filterLeft="",filterRight=""):
        f= open("html/runpythonbiogeme.html")
        page = f.read()

        modelFiles = self.getDragableTable(".py",filterLeft,sort, inspect.stack()[0][3])
        dataFiles = self.getDragableTable(".dat",filterRight,sort, inspect.stack()[0][3])
        navBar = self.getNavBar(2,1)
        footer = self.getFooter()

        page = page.replace("MODELS", modelFiles)
        page = page.replace("DATA", dataFiles)
        page = page.replace("DYNAMIC_NAV_BAR", navBar)
        page = page.replace("DYNAMIC_FOOTER", footer)
        
        return page
    runpb.exposed = True

    #page to choose the model and data to launch with biogeme
    def runb(self, sort=DEFAULT_SORT,filterLeft="",filterRight=""):
        f= open("html/runbiogeme.html")
        page = f.read()
        
        modelFiles = self.getDragableTable(".mod",filterLeft,sort, inspect.stack()[0][3])
        dataFiles = self.getDragableTable(".dat",filterRight,sort, inspect.stack()[0][3])
        navBar = self.getNavBar(2,2)
        footer = self.getFooter()

        page = page.replace("MODELS", modelFiles)
        page = page.replace("DATA", dataFiles)
        page = page.replace("DYNAMIC_NAV_BAR", navBar)
        page = page.replace("DYNAMIC_FOOTER", footer)

        return page
    runb.exposed = True

    #page to choose the model and data to launch with biosim
    def runbs(self, sort=DEFAULT_SORT,filterLeft="",filterRight=""):
        f = open("html/runbiosim.html")
        page = f.read()
        
        modelFiles = self.getDragableTable(".mod",filterLeft,sort, inspect.stack()[0][3])
        dataFiles = self.getDragableTable(".dat",filterRight,sort, inspect.stack()[0][3])
        navBar = self.getNavBar(2,3)
        footer = self.getFooter()

        page = page.replace("MODELS", modelFiles)
        page = page.replace("DATA", dataFiles)
        page = page.replace("DYNAMIC_FOOTER", footer)
        page = page.replace("DYNAMIC_NAV_BAR", navBar)
        
        return page
    runbs.exposed = True

    #page that displays the console output of the command
    def console(self, cmd="0", mod="", dat=""):
        started = self.callBiogeme(cmd, mod, dat)
        f = open("html/console.html")
        page = f.read()
        if not self.busy and self.lastModelUsed != "":
            page = page.replace("RESULTS",self.linkToResults())
        else:
            page = page.replace("RESULTS", "")
        

        m = ""
        if started:
        
            if mod != "" and dat != "":
                m = self.prepareAlert("Started with [" + mod + "] and [" + dat+ "]", "success")
            elif cmd.startswith("ping"):
                m = self.prepareAlert("ping self "+ cmd[-2:] +" times for test prupose", "info")
            else:
                pass
                #m = self.prepareAlert("started one or two empty", "info")
        else:
            if mod != "" and dat != "" and (mod != self.lastModelUsed or dat != self.lastDataUsed):
                m = self.prepareAlert(self.lastModelUsed + " " + self.lastDataUsed + " already running", "danger")
            elif mod != "" and dat != "" and mod == self.lastModelUsed and dat == self.lastDataUsed:
                m = self.prepareAlert("[" + self.lastModelUsed + "] and [" + self.lastDataUsed + "] already running", "warning")
            elif "ping" == self.lastModelUsed and "ping" == self.lastDataUsed:
                m = self.prepareAlert("ping already running", "warning")
            else:
                m = self.prepareAlert("Could not start with ["+mod+"] and ["+dat+"], please verify your files.", "danger")

        page = page.replace("MESSAGE", m)

        footer = self.getFooter()
        page = page.replace("DYNAMIC_FOOTER", footer)
        navBar = self.getNavBar(3)
        page = page.replace("DYNAMIC_NAV_BAR", navBar)

        return page
    console.exposed = True


    #page to display the html and latex results
    #can be filtered and/or sorted
    def results(self,filter="",sort=DEFAULT_SORT):
        f = open("html/results.html")
        page = f.read()


        page = page.replace("HTML_RESULTS", self.getDragableTable(".html",filter,sort, inspect.stack()[0][3]) )
#        page = page.replace("TEX_RESULTS", self.getDragableTable(".tex",filter,sort, inspect.stack()[0][3]) )
        page = page.replace("LOG_RESULTS", self.getDragableTable(".log",filter,sort, inspect.stack()[0][3]) )

        footer = self.getFooter()
        page = page.replace("DYNAMIC_FOOTER", footer)
        navBar = self.getNavBar(4)
        page = page.replace("DYNAMIC_NAV_BAR", navBar)
         
        return page
    results.exposed = True

    #page to display the documentation
    def doc(self):
        f = open("html/doc.html")
        page = f.read()
        
        if os.path.isfile("html/doxygen/html/annotated.html"):
            page = page.replace("DYNAMIC_DOC_PATH", "doxygen/annotated.html")
        else:
            page = page.replace("DYNAMIC_DOC_PATH", "about:blank")

        footer = self.getFooter()
        page = page.replace("DYNAMIC_FOOTER", footer)
        navBar = self.getNavBar(5)
        page = page.replace("DYNAMIC_NAV_BAR", navBar)
        return page

    doc.exposed = True

    #display a quick tutorial in a bootstrap carousel.
    #the parts about the settings page are hidden if the settings are hidden
    def help(self):
        f = open("html/help.html")
        page = f.read()
        
        footer = self.getFooter()
        page = page.replace("DYNAMIC_FOOTER", footer)
        navBar = self.getNavBar(6)
        page = page.replace("DYNAMIC_NAV_BAR", navBar)

        if not SHUTDOWNABLE:
            re.DOTALL
            a = re.sub(r"(<!--SHUTDOWNABLE_BEGIN-->).*?(<!--SHUTDOWNABLE_END-->)", "", page, 0,re.DOTALL)
            print(a)
            return a

        return page
    help.exposed = True

    #settings page, where we can set the :
    # sort: default sorting when some files are listed.
    # data: working directory, can be absolute or relative
    # port: choose the port on which the webui runs
    # visi: 0.0.0.0 or 127.0.0.1, the first allows every ip address to access the webui, the second only from local ip
    def settings(self,sort='',data='',port='',visi=''):
        global DATA_FOLDER
        f = open("html/settings.html")
        page = f.read()


        configBiogemeChanged = False
        if sort != '' and sort != config["biogemeWebUI"]["defaultSort"]:
            config.set("biogemeWebUI","defaultSort",sort)
            DEFAULT_SORT = sort
            configBiogemeChanged = True

        if data != '' and data != config["biogemeWebUI"]["datafolder"]:
            config.set("biogemeWebUI","datafolder",data)
            DATA_FOLDER = os.path.join(os.getcwd(), data)
            print("********************** New data folder: " + DATA_FOLDER + " - " +os.getcwd())
            if not os.path.isdir(DATA_FOLDER):
                print("- Making directory")
                os.makedirs(DATA_FOLDER)
            configBiogemeChanged = True


        configServerChanged = False
        if visi != '' and visi != configServer["global"]["server.socket_host"][1:-1]:
            print("**********************" + visi + " " + configServer["global"]["server.socket_host"])
            configServer.set("global", "server.socket_host", '"' + str(visi) + '"')
            configServerChanged = True

        if port != '' and port != configServer["global"]["server.socket_port"]:
            print("**********************" + port + " " + configServer["global"]["server.socket_port"])
            configServer.set("global", "server.socket_port", str(port))
            configServerChanged = True


        page = page.replace("DYN_"+config["biogemeWebUI"]["defaultsort"],"selected")
        page = page.replace("DYN_DATA_FOLDER",config["biogemeWebUI"]["datafolder"])

        page = page.replace("DYN_"+configServer["global"]["server.socket_host"][1:-1],"selected")
        page = page.replace("DYN_SERVER_PORT", configServer["global"]["server.socket_port"])



        msg = ""
        if configServerChanged or configBiogemeChanged:
            msg = self.prepareAlert("Settings updated.")

        if port != '' and port != configServer["global"]["server.socket_port"]:
            msg = self.prepareAlert("You changed the port to " + port + ", please correct the url.", "warning")

        if onWindows and configServerChanged:
            msg += self.prepareAlert("New server settings will take effect on next start.","danger")


        page = page.replace("DYN_MESSAGE",msg)
        footer = self.getFooter()
        page = page.replace("DYNAMIC_FOOTER", footer)
        navBar = self.getNavBar(7)
        page = page.replace("DYNAMIC_NAV_BAR", navBar)


        #check if restart is needed
        if configBiogemeChanged:
            config.write(open("biogemeWebUI.conf", "w"))
            if not onWindows:
                print("restart")
                cherrypy.engine.restart()

        if configServerChanged:
            configServer.write(open("server.conf", "w"))
            if not onWindows:
                cherrypy.engine.restart()

        
        return page
    settings.exposed = SHUTDOWNABLE

    

    #page to shutdown the webui, will display a progress bar that is totally random,
    # but it shows that the shutdown progress takes some time.
    def shutdown(self):
        f = open("html/closing.html")
        page = f.read()

        #alert = self.prepareAlert("Closing","info")
        #page = page.replace("DYNAMIC_MAIN", alert)

        navBar = self.getNavBar(8)
        page = page.replace("DYNAMIC_NAV_BAR", navBar)
        footer = self.getFooter()
        page = page.replace("DYNAMIC_FOOTER", footer)

        cherrypy.engine.exit()
        return page
    shutdown.exposed = SHUTDOWNABLE


    def version(self):
        rep = b'biogemeWebUI<br />'
        v = subprocess.Popen(["pythonbiogeme", "--version"], stdout=subprocess.PIPE,shell=onWindows)
        rep += v.stdout.readline()
        return rep
    version.exposed = True

    #========================================================================

    




if __name__ == '__main__':
    # CherryPy always starts with app.root when trying to map request URIs
    # to objects, so we need to mount a request handler root. A request
    # to '/' will be mapped to biogemeWebUI().index().
    confdict = {'/bootstrap':    {'tools.staticdir.on': True, 
                        'tools.staticdir.dir': BS_FOLDER},
            '/jasny-bootstrap':{'tools.staticdir.on': True, 
                        'tools.staticdir.dir': JBS_FOLDER},
            '/data':    {'tools.staticdir.on': True, 
                        'tools.staticdir.dir': DATA_FOLDER},
            '/img':    {'tools.staticdir.on': True, 
                        'tools.staticdir.dir': IMG_FOLDER},
            '/css':     {'tools.staticdir.on': True,
                        'tools.staticdir.dir': CSS_FOLDER},
            '/doxygen':     {'tools.staticdir.on': True,
                        'tools.staticdir.dir': DOXY_FOLDER},
            '/favicon.ico':{'tools.staticfile.on': True,
                          'tools.staticfile.filename': ICONE}
            }


    app = cherrypy.tree.mount(BiogemeWebUI(), '/', confdict)

    cherrypy.config.update("server.conf")



    if hasattr(cherrypy.engine, "signal_handler"):
        cherrypy.engine.signal_handler.subscribe()
    if hasattr(cherrypy.engine, "console_control_handler"):
        cherrypy.engine.console_control_handler.subscribe()

    
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.settimeout(0.2)
        s.connect(("localhost", int(configServer["global"]["server.socket_port"])))
    except socket.timeout:
        #windows timeout when port is free
        cherrypy.engine.start()
    except ConnectionRefusedError:
        #unix get a connection refused if port is free
        cherrypy.engine.start()
    else:
        print("Port " + configServer["global"]["server.socket_port"] + " already in use.")
        quit()

    finally:
        if len(sys.argv) > 1:
            command = str(sys.argv[1])
            if command == "openBrowser":
                print("Open browser on port "+ configServer["global"]["server.socket_port"])
                webbrowser.open("http://127.0.0.1:" + configServer["global"]["server.socket_port"] +"/index")
            else:
                print("Unknow command: " + command)


    cherrypy.engine.block()
    print("BiogemeWebUI closed.")
