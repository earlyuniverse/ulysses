#!/usr/bin/env python3
import glob
import re

# first identify all AUTHORS files
AuthorFiles = glob.glob("*/AUTHORS")

print ("Identified {} contribs with AUTHORS files\n".format(len(AuthorFiles)))

# create an initial empty set for emails
emails = set()

# populate it
for file in AuthorFiles:
    with open(file,'r') as f:
        for line in f:
            line = line.strip()
            # pull out all email addresses
            element = r'[^@ <>(),:]+'
            m = re.findall('('+element+'@'+element+')',line)
            for email in m: 
                email = re.sub(r'\.$','',email)
                emails.add(email)

# print them in a format that we can easily check
print("\n".join(emails),"\n")

# print them in a format that works for a mail program...
print(",".join(emails))
