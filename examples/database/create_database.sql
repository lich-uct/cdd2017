create extension if not exists rdkit;

drop schema if exists drugbank cascade;
create schema drugbank;

-- main table
create table drugbank (
       id serial primary key,
       smiles text,
       molinstance mol
);
