#! /usr/bin/env python
print("TEST!")
import socket
import sys

print "Loading popcoen modules..."

import popcoenMain

print "Done."

import tempfile
from cStringIO import StringIO
import traceback


# Create a TCP/IP socket
sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

# Bind the socket to the port
# start server locally:
server_address = ('localhost', 10888)

# through WWW works similarly. Firewall needs to allow it though:
# servername = socket.gethostname()
# server_address = (servername, 10888)


print 'Starting up Popcoen-server on %s port %s' % server_address
sock.bind(server_address)



# Listen for incoming connections
sock.listen(1)

while True:
    # Wait for a connection
    print '\n-------------------\nPopcoen-server is listening, i.e. waiting for incoming connection...'
    connection, client_address = sock.accept()

    try:
        print 'Connection from', client_address

        # Receive the data in small chunks and retransmit it
	TOTDATA = ""
        while True:
            data = connection.recv(1000)
            TOTDATA += data

            #print >>sys.stderr, 'received "%s"' % data
            if TOTDATA.find("SENDEND") != -1:
        	if (TOTDATA[0:8] != "PDBstart") or (TOTDATA[-13:] != "PDBendSENDEND"):
        	    print "Error Transmission from client to server"
        	    BACK = "ERROR transmission_RESPONSEEND"
        	else:

		    print 'No more data from', client_address
		    print "Number bytes received:", len(TOTDATA)

        	    pdbcontent = TOTDATA[8:-13]

        	    t = tempfile.NamedTemporaryFile(mode='w+b', bufsize=-1, suffix='.pdb', prefix='tmpfile_mdtrajInput_', dir="./")
        	    t.write(pdbcontent)
        	    t.seek(0)

        	    pdbfilename = t.name
        	    print "Data saved to tmpfile",pdbfilename,"and will be given to popcoen calculation"



        	    #We just redirect stdout so that Main popcoen program does not have to be changes
                    old_stdout = sys.stdout
                    sys.stdout = mystdout = StringIO()

                    try:
			#calculate popcoen entropies here
        	    	popcoenMain.Popcoen_all___from_structure_to_prediction(pdbfilename)
                        BACK = mystdout.getvalue()
                        BACK = "RESPONSESTART%sRESPONSEEND" %BACK
                    except:
                    	errormessage = traceback.format_exc()  # catch error message
                        BACK = "RESPONSESTARTERROR occured:\n%sRESPONSEEND" %errormessage

                    sys.stdout = old_stdout
                    print "Sending back:\n",BACK[0:200],"[...]",BACK[-100:]



		# send back data or error
		connection.sendall(BACK)
		break
            
    finally:
        # Clean up the connection
        connection.close()
