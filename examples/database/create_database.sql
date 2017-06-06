create extension if not exists rdkit;

drop schema if exists drugbank cascade;
create schema drugbank;

-- main table
create table drugbank.drugbank (
       id serial primary key,
       smiles text,
       molinstance mol,
       morganfp bfp,
       hba integer,
       hbd integer,
       mw float,
       logp float
);

-- drug groups table
create table drugbank.druggroup (
       id serial primary key,
       groupname text
);

-- drug group <--> drug many2many connection table
create table drugbank.druggroup2drugbank (
       drugbank_id integer references drugbank.drugbank(id),
       druggroup_id integer references drugbank.druggroup(id)
);
