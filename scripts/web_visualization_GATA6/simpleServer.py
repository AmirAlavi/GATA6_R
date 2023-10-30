# This script only works for python 3.x
 
from http.server import HTTPServer
from http.server import CGIHTTPRequestHandler as handler
import cgitb; cgitb.enable()  ## This line enables CGI error reporting
import sys
    
server_address = ("", int(sys.argv[1]))
handler.cgi_directories = ["/ht-bin"]
 
httpd = HTTPServer(server_address, handler)
httpd.serve_forever()