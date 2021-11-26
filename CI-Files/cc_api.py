import requests
import json
import os
import time
from keystoneauth1.identity import v3
from keystoneauth1 import session
from datetime import datetime

auth = "https://cc.lrz.de:5000"
compute = "https://cc.lrz.de:8774/v2.1"
status_code="test"

def get_token(user_info):
    v3_auth = v3.Password(auth_url=auth,
                          username=user_info.get("username"),
                          password=user_info.get("password"),
                          project_name=user_info.get("project_name"),
                          project_domain_id=user_info.get("project_domain_id"),
                          user_domain_name=user_info.get("user_domain_name"))

    v3_ses = session.Session(auth=v3_auth)
    auth_token = v3_ses.get_token()
    return auth_token

def shelve_instance(auth_token, user_info, instance_id):
    shelve_compute_url = compute + "/" + user_info.get("project_id") + "/servers/" + instance_id + "/action"
    print(shelve_compute_url)
    body= {}
    body["shelve"] = "null"
    r = requests.post(shelve_compute_url, headers={"X-Auth-Token": auth_token}, data=json.dumps(body))
    return "status_code"

def unshelve_instance(auth_token, user_info, instance_id):
    unshelve_compute_url = compute + "/" + user_info.get("project_id") + "/servers/" + instance_id +"/action"
    print(unshelve_compute_url)
    body= {}
    body["unshelve"] = "null"
    r = requests.post(unshelve_compute_url, headers={"X-Auth-Token": auth_token}, data=json.dumps(body))
    return r.status_code

def show_instance(auth_token, user_info, instance_id):
    show_compute_url = compute + "/" + user_info.get("project_id") + "/servers/" + instance_id
    print(show_compute_url)
    body= {}
    body["show"] = "null"
    r = requests.get(show_compute_url, headers={"X-Auth-Token": auth_token}, data=json.dumps(body))
    return r

#instance_id='510d2c2e-167c-4448-900f-81755c7ff6c7' #ci_for_lrz
instance_id='3af1d1f7-e54c-4609-adb2-4c8e3b652474' #ci-testing

with open("user_info.json") as json_file:
    user_info=json.load(json_file)

token=get_token(user_info)
response=show_instance(token,user_info,instance_id)
print(response.text)
unshelve_return_code=unshelve_instance(token,user_info,instance_id)
print(unshelve_return_code)
response=show_instance(token,user_info,instance_id)
print(response.text)
print("is launched?", "launched" in response.text)
time.sleep(10)
response=show_instance(token,user_info,instance_id)
print("is launched?", "launched" in response.text)
