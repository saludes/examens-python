#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Filename:   enviar-examens.py
Author:     Rafel Amer (rafel.amer AT upc.edu)
Copyright:  Rafel Amer 2018
Disclaimer: This code is presented "as is" and it has been written to
            generate random models of exams for the subject of Linear
            Algebra at ESEIAAT, Technic University of Catalonia
License:    This program is free software: you can redistribute it and/or modify
            it under the terms of the GNU General Public License as published by
            the Free Software Foundation, either version 3 of the License, or
            (at your option) any later version.

 	        See https://www.gnu.org/licenses/
"""

SCOPES = ['https://www.googleapis.com/auth/gmail.readonly','https://www.googleapis.com/auth/gmail.send']

import pickle
import mimetypes
import os.path
import base64
import mimetypes
import sys
import re
import unidecode
from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders
from googleapiclient.errors import HttpError
from optparse import OptionParser

parser = OptionParser()
parser.add_option("--estudiants",dest="estudiants",default=None)
parser.add_option("--subject",dest="subject",default=None)
parser.add_option("--message",dest="message",default=None)
parser.add_option("--sender",dest="sender",default=None)
parser.add_option("--ajuda",action="store_true",dest="ajuda",default=False)
parser.add_option("--solucions",action="store_true",dest="solucions",default=False)
(options,args) = parser.parse_args()

def create_message(sender,to,subject,message_text,files):
    """Create a message for an email.
    Args:
        sender: Email address of the sender.
        to: Email address of the receiver.
        subject: The subject of the email message.
        message_text: The text of the email message.
        files: List of the paths to the files to be attached.
    Returns:
    An object containing a base64url encoded email object.
    """
    message = MIMEMultipart()
    message['to'] = to
    message['from'] = sender
    message['subject'] = subject

    msg = MIMEText(message_text)
    message.attach(msg)

    for file in files:
        content_type, encoding = mimetypes.guess_type(file)
        if content_type is None or encoding is not None:
            content_type = 'application/octet-stream'
        main_type, sub_type = content_type.split('/', 1)
        fp = open(file, 'rb')
        msg = MIMEBase(main_type, sub_type)
        msg.set_payload(fp.read())
        fp.close()
        encoders.encode_base64(msg)
        filename = os.path.basename(file)
        msg.add_header('Content-Disposition', 'attachment', filename=filename)
        message.attach(msg)
    return {'raw': base64.urlsafe_b64encode(bytes(message.as_string().encode('utf-8'))).decode("ascii")}

def send_message(service, user_id, message):
    """
    Send an email message.
    Args:
        service: Authorized Gmail API service instance.
        user_id: User's email address. The special value "me"
            can be used to indicate the authenticated user.
        message: Message to be sent.

    Returns:
        Sent Message.
    """
    try:
        message = service.users().messages().send(userId=user_id, body=message).execute()
        print ('Message Id: %s' % message['id'])
    except HttpError as error:
        print ('An error occurred: %s' % error)

parser = OptionParser()
parser.add_option("--estudiants",dest="estudiants",default=None)
parser.add_option("--subject",dest="subject",default=None)
parser.add_option("--message",dest="message",default=None)
parser.add_option("--sender",dest="sender",default=None)
parser.add_option("--ajuda",action="store_true",dest="ajuda",default=False)
parser.add_option("--solucions",action="store_true",dest="solucions",default=False)
(options,args) = parser.parse_args()

HOME = os.path.expanduser('~')
est = options.estudiants
fitxer = options.message
regex = re.compile('^\s*#.$',re.IGNORECASE)
estudiants = []
sender, subject = options.sender,options.subject
if sender is None or subject is None:
    print("S'ha d'especificar l'assumpte i l'emisor dels correus")
    sys.exit(0)

try:
    with open(est) as f:
        for line in f:
            line = line.rstrip()
            if regex.match(line):
                continue
            try:
                data = line.split(':')
                estudiants.append({'nom' : data[0],'cognoms' : data[1],'email' : data[2]})
            except:
                continue
    f.close()
except:
    print("Error de lectura del fitxer d'estudiants")
    sys.exit(0)

try:
    with open(fitxer,'r') as f:
        message = f.read()
except:
    print("Error de lectura del fitxer amb el missatge del correu")
    sys.exit(0)

try:
    with open(f"{HOME}/credentials/token.pickle", 'rb') as token:
        creds = pickle.load(token)
except:
    print(f"Error en llegir el fitxer {HOME}/credentials/token.pickle")
    sys.exit(0)

if not creds or not creds.valid:
    if creds and creds.expired and creds.refresh_token:
        creds.refresh(Request())
    else:
        flow = InstalledAppFlow.from_client_secrets_file(f"{HOME}/credentials/credentials.json", SCOPES)
        creds = flow.run_local_server(port=0)
        with open(f"{HOME}/credentials/token.pickle", 'wb') as token:
            pickle.dump(creds, token)

service = build('gmail', 'v1', credentials=creds)

for e in estudiants:
    relacio = relacio = {'COGNOMS' : e['cognoms'], 'NOM' : e['nom']}
    m = message
    for k,v in relacio.items():
        m = m.replace(k,v)
    filename = f"tex/{e['cognoms']}-{e['nom']}".lower().replace(' ','-')
    filename = unidecode.unidecode(filename)
    if options.solucions:
        filename += "-solucio.pdf"
    else:
        filename += ".pdf"
    correu = create_message(sender,e['email'],subject,m,[filename])
    send_message(service,'me', correu)
