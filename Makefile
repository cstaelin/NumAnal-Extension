# If this Makefile is not being run from the NetLogo-6.x.x/app/extensions/numanal
# directory, NETLOGO must be set before the make to the location of the NetLogo
# package.  E.g., in Cygwin this would be something like 
# export NETLOGO = "c:/Program Files/NetLogo 6.1.1"

ifeq ($(origin JAVA_HOME), undefined)
	JAVA_HOME = /usr
endif

ifneq (,$(findstring CYGWIN,$(shell uname -s)))
  COLON=\;
else
  COLON=:
endif

ifeq ($(origin NETLOGO), undefined)
	NETLOGO = ../../..
endif

# NetLogo.jar files are now modified by the NetLogo version number,
# e.g., netlogo-6.1.0.jar. Thus we shell to the "find" command which 
# looks in the app folder to get the .jar file name. I tried to 
# accomplish this with the wildcard function, but it doesn't seem to 
# handle cases where there are blanks in file/directory names, as is 
# common in Windows.
NETLOGO_JAR := "$(shell find "$(NETLOGO)"/app -name netlogo-*.jar)"
#NETLOGO_JAR := $(wildcard $(NETLOGO)/app/netlogo-*.jar)
JAVAC := "$(JAVA_HOME)/bin/javac"
JAVAJAR := "$(JAVA_HOME)/bin/jar"
SRCS := $(wildcard src/*.java)

NumAnalExtension.zip: numanal.jar Jama-1.0.3.jar commons-math3-3.6.1.jar PalMathLibrary.jar README.md license.md Makefile src manifest.txt NumAnal-v3.3.0.pdf Examples
	rm -rf numanal
	mkdir numanal
	cp -rp numanal.jar Jama-1.0.3.jar commons-math3-3.6.1.jar PalMathLibrary.jar README.md license.md Makefile src manifest.txt NumAnal-v3.3.0.pdf Examples numanal
	zip -rv NumAnalExtension.zip numanal
	rm -rf numanal

numanal.jar: $(SRCS) Jama-1.0.3.jar commons-math3-3.6.1.jar PalMathLibrary.jar Makefile manifest.txt
	rm -rf classes
	mkdir -p classes
	$(JAVAC) -g -deprecation -Xlint:all -Xlint:-serial -Xlint:-path -encoding us-ascii -source 1.8 -target 1.8 -classpath $(NETLOGO_JAR)$(COLON)Jama-1.0.3.jar$(COLON)commons-math3-3.6.1.jar$(COLON)PalMathLibrary.jar -d classes $(SRCS)
	$(JAVAJAR) cmf manifest.txt numanal.jar -C classes .
	rm -rf classes


