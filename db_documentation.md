#something

create table address
(
    address_id serial
        primary key,
    country    varchar(255),
    region     varchar(255),
    city       varchar(255),
    district   varchar(255),
    post_code  varchar(255),
    street     varchar(255),
    building   varchar(255),
    room       varchar(255)
);
