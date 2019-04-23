import os
os.system("auth-tokens-refresh")
for i in range(10):
    os.system("abaco deploy")
