import json

config = json.load(open('config.json'))

print("LR is " + config["mlp"]["lr"]);
