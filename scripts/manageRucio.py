#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import os
import argparse
import yaml
import subprocess
from rucio.client import Client

PREFIX = "/pnfs/psi.ch/cms/trivcat"

def yesno(question):
    resp = input(question + " (yes/no) ")
    return (resp.lower() == "yes") or (resp.lower() == "y")

def create_container(container):
    subprocess.run(["rucio", "add-container", container])

def create_rule(container, site, comment):
    subprocess.run(["rucio", "add-rule", "--ask-approval",
                    "--comment", comment, container, "1", site])

def list_rules(user):
    # subprocess.run(["rucio", "list-rules", "--account", user])
    client = Client()
    resp = client.list_account_rules(user)
    for rule in resp:
        print(rule)

def rule_status(rule):
    subprocess.run(["rucio", "rule-info", rule])
    # client = Client()
    # resp = client.get_replication_rule(rule)

def get_container_content(container):
    scope,name = container.split(":")
    client = Client()
    resp = client.list_content(scope, name)
    return list(r["name"] for r in resp)

def list_container(container):
    # subprocess.run(["rucio", "list-content", container])
    for smp in get_container_content(container):
        print(smp)

def check_files(container):
    scope,name = container.split(":")
    client = Client()
    resp = client.list_files(scope, name)
    import os.path
    print("The following files are missing from the container:")
    for r in resp:
        full_path = f"{PREFIX}/{r['name']}"
        if not os.path.isfile(full_path):
            print(full_path)


def add_to_container(container, sample_list):
    did_list = [ "cms:" + smp for smp in sample_list ]
    result = subprocess.run(["rucio", "attach", container] + did_list)
    return result.returncode == 0

def remove_from_container(container, sample_list):
    did_list = [ "cms:" + smp for smp in sample_list ]
    result = subprocess.run(["rucio", "detach", container] + did_list)
    return result.returncode == 0

def check_sample_exists(sample):
    # For some reason this sometimes fails:
    # client = Client()
    # resp = client.list_replicas([{"scope": "cms", "name": sample}])
    # try:
    #     next(resp)
    # except StopIteration:
    #     return False
    # else:
    #     print(f"Sample {sample} exists!")
    #     return True
    result = subprocess.run(["rucio", "list-dataset-replicas", "cms:" + sample],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return result.returncode == 0

def sync_container(container, sample_template, eras=None):
    sample_list_from_yaml = set()
    for sample_name,sample_config in sample_template.items():
        if "db" in sample_config:
            if eras is not None and sample_config["era"] not in eras:
                continue
            sample_list_from_yaml.add(sample_config["db"].split("das:")[1])
        elif "dbs" in sample_config:
            for era,das_key in sample_config["dbs"].items():
                if eras is not None and era not in eras:
                    continue
                sample_list_from_yaml.add(das_key.split("das:")[1])

    sample_list_from_container = set(get_container_content(container))

    to_add = sample_list_from_yaml.difference(sample_list_from_container)
    to_remove = sample_list_from_container.difference(sample_list_from_yaml)

    if to_add:
        for smp in to_add:
            if not check_sample_exists(smp):
                raise RuntimeError(f"ERROR: sample {smp} does not exist!")
            print(smp)
        if yesno(f"\nThe above list of samples is going to be ADDED to the container {container}, would you like to continue?"):
            success = add_to_container(container, to_add)
            if success:
                print("\nSuccess!")
            else:
                print("\nSomething went wrong!")
        else:
            print("\nAbort!")
        print("\n")

    if to_remove:
        for smp in to_remove:
            print(to_remove)
        if yesno(f"\nThe above list of samples is going to be REMOVED from the container {container}, would you like to continue?"):
            success = remove_from_container(container, to_remove)
            if success:
                print("\nSuccess!\n\n")
            else:
                print("\nSomething went wrong!\n\n")
        else:
            print("\nAbort!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Utilities to manage Rucio rules and containers')
    parser.add_argument('--samples', type=str, nargs="*", help='Sample template YAML file')
    parser.add_argument('--container', type=str, help='Specify container name (with user.<user> prefix)')
    parser.add_argument('--rule', type=str, help='Specify rule hash')
    parser.add_argument('--eras', type=list, help='Filter eras in sample list (careful about removal from container!!)')
    parser.add_argument('--user', type=str, help='Specify username (by default take current user)', default=os.getenv("USER"))
    parser.add_argument('--site', type=str, help='Grid site for replica storage')
    parser.add_argument('--comment', type=str, help='Comment for rule approval')
    parser.add_argument('--create-container', action='store_true', help='Create container')
    parser.add_argument('--sync', action='store_true', help='Sync container with sample template list')
    parser.add_argument('--create-rule', action='store_true', help='Create rule for given container at given site')
    parser.add_argument('--list-rules', action='store_true', help='List rules associated with a given account')
    parser.add_argument('--rule-status', action='store_true', help='Show status of given rule')
    parser.add_argument('--list-container', action='store_true', help='Show content of given container')
    parser.add_argument('--check-files', action='store_true', help='Check if all files in the container are present locally')
    options = parser.parse_args()

    if options.create_container:
        assert(options.container is not None)
        create_container(options.container)
    if options.create_rule:
        assert(options.container is not None)
        assert(options.site is not None)
        assert(options.comment is not None)
        create_rule(options.container, options.site, options.comment)
    if options.list_rules:
        list_rules(options.user)
    if options.rule_status:
        assert(options.rule is not None)
        rule_status(options.rule)
    if options.list_container:
        assert(options.container is not None)
        list_container(options.container)
    if options.sync:
        assert(options.container is not None)
        assert(options.samples is not None)
        sample_template = dict()
        for path in options.samples:
            with open(path) as fh:
                sample_template.update(yaml.safe_load(fh))
        sync_container(options.container, sample_template, options.eras)
    if options.check_files:
        assert(options.container is not None)
        check_files(options.container)
