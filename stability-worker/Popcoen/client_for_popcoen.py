#! /usr/bin/env python

import socket
import sys

from timeit import default_timer as timer


pdbContent = sys.stdin.read()


start = timer()

# Create a TCP/IP socket
sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

# Connect the socket to the port where the server is listening
# locally:
server_address = ('localhost', 10888)
# globally would works similarly:
#server_address = ('Obelix.ffn.ub.es', 10888)

print 'INFO: Connecting to popcoen (at %s port %s)' % server_address
sock.connect(server_address)

try:
    
    #f = open("1q2w.pdb")
    message = "PDBstart" + pdbContent + "PDBendSENDEND"
    #f.close()
    # Send data
    #message = 'This is the message.  It will be repeated.'
    print 'INFO: Sending pdb-file to popcoen server.'
    sock.sendall(message)

    # Look for the response
    GOTBACK = ""
    while GOTBACK.find("RESPONSEEND") == -1:
        data = sock.recv(1000)
	GOTBACK += data

    if (GOTBACK[0:13] != "RESPONSESTART") or (GOTBACK[-11:] != "RESPONSEEND"):
	print "Transmission Error occured"
    else:
	Result = GOTBACK[13:-11]
	#print "Result obtained:\n>%s<" % Result
	print Result

finally:
    print 'INFO: closing socket to popcoen server.'
    sock.close()

    end = timer()
    print "INFO: Elapsed time for calculation and socket-communication= %.1lf s" % (end-start)
