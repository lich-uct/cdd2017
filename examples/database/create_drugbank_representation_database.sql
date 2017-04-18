create extension if not exists rdkit;

drop schema if exists drugbank cascade;
create schema drugbank;

-- main table
create table drugbank.drugbank (
       id text primary key,
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
       drugbank_id text references drugbank.drugbank(id),
       druggroup_id integer references drugbank.druggroup(id)
);

-- product name --> API one to many
create table drugbank.products (
       productname text, -- can't really be unique here
       drugbank_id text references drugbank.drugbank(id)
);
