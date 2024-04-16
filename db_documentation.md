## 1. Таблица address

```sql 
create table address -- всё об адресах
( 
    address_id serial 
        primary key, -- айдишник, первичный ключ 
    country    varchar(255), -- страна 
    region     varchar(255), -- регион или область 
    city       varchar(255), -- город 
    district   varchar(255), -- район 
    post_code  varchar(255), -- почтовый индекс 
    street     varchar(255), -- улица 
    building   varchar(255), -- номер дома 
    room       varchar(255) -- номер офиса
);

alter table address
    owner to hrtools;

grant select on address to grafana;
```

## 2. Таблица attachable_permission
``` sql
create table attachable_permission -- должности, сокращения и т.п.
( 
    attachable_permission_id serial 
        primary key, -- айдишник, первичный ключ
    name                     varchar(255), -- имя в системе по типу attachable.god
    description              varchar(255), -- описание имени, например, Партнёр
    shortcut                 varchar(255), -- сокращение (пока они такие же как и имена) 
    "order"                  integer -- порядок расположения на сайте
); 

alter table attachable_permission 
    owner to hrtools; 

grant select on attachable_permission to grafana; 
```

## 3. Таблица attachable_permission_grant
```sql
create table attachable_permission_grant 
( 
    attachable_permission_grant_id serial 
        primary key, -- айдишник, первичный ключ
    auth_item_name                 varchar(255) not null, --  имя в системе по типу attachable.god, аналогично attachable_permission
    attachable_permission_id       integer      not null 
        references attachable_permission -- ссылается на айдишник (attachable_permission_id) в табличке attachable_permission
            on update cascade on delete cascade 
); 
alter table attachable_permission_grant
    owner to hrtools; 

create index attachable_permission_grant_attachable_permission_id_idx
    on attachable_permission_grant (attachable_permission_id);

grant select on attachable_permission_grant to grafana; 
```

## 4. Таблица table auth_assignment
``` sql
-- auto-generated definition 
create table auth_assignment -- идентификация пользователей на платформе(?) auth_item
( 
    item_name  varchar(64) not null 
        references auth_item -- ссылается на статус (name) в таблице 
            on update cascade on delete cascade, -- статус, например employee, registered; ссылается статус (name) в табличке auth_item 
    user_id    integer     not null, -- айдишник 
    created_at integer, -- когда был создан (формат не временной?) 
    primary key (item_name, user_id) 
);
```

## 5. Таблица auth_item
```sql
-- auto-generated definition 
create table auth_item 
( 
    name        varchar(64) not null 
        primary key, -- статус по типу site.access, /, /site/*, rbac.edit
    type        smallint    not null, -- тип, пока их всего 2
    description text, -- описание каждого статуса, что и кому доступно
    rule_name   varchar(64) 
                            references auth_rule 
                                on update cascade on delete set null,
-- название правил по типу OwnVacancyRule ?,
-- ссылается на табличку auth_rule, где все правила описаны
    data        bytea, -- видимо, какие-то данные, но пока всё пустое ?
    created_at  integer, -- когда создавалось (формат числовой, не даты)
    updated_at  integer -- когда обновлено (формат числовой, не даты)
); 

alter table auth_item 
    owner to hrtools; 

create index "idx-auth_item-type" 
    on auth_item (type); 

grant select on auth_item to grafana; 

grant select on auth_item to clickhouse; 


alter table auth_assignment 
    owner to hrtools; 

create index "idx-auth_assignment-user_id" 
    on auth_assignment (user_id); 

grant select on auth_assignment to grafana; 

grant select on auth_assignment to clickhouse; 

```

## 6. Таблица auth_item_child

``` sql
-- У кого какие дети-процессы?
-- auto-generated definition
create table auth_item_child
(
    parent varchar(64) not null
        references auth_item
            on update cascade on delete cascade, -- ссылается на статусы (name)  в auth_item
    child  varchar(64) not null
        references auth_item
            on update cascade on delete cascade, -- ссылается на статусы (name)  в auth_item
    primary key (parent, child) -- составной первичный ключ
);

alter table auth_item_child
    owner to hrtools;

grant select on auth_item_child to grafana;

grant select on auth_item_child to clickhouse;

```

## 7. Таблица auth_rule

``` sql
-- auto-generated definition
create table auth_rule
(
    name       varchar(64) not null 
        primary key, -- название правил
    data       bytea, -- текстовый данные, видимо словарь ?
    created_at integer, -- когда создавалось (формат числовой, не даты)
    updated_at integer -- когда обновлялось (формат числовой, не даты)
);

alter table auth_rule
    owner to hrtools;

grant select on auth_rule to grafana;

grant select on auth_rule to clickhouse;

```

## 8. Таблица candidate

``` sql
-- auto-generated definition
create table candidate -- сборная солянка по кандидатам
(
    candidate_id        serial
        primary key, -- первичный ключ, айдишник кандидата
    first_name          varchar(255), -- имя
    last_name           varchar(255), -- фамилия
    middle_name         varchar(255), -- отчество
    birth_year          integer, -- год рождения
    position            varchar(255), -- позиция, на которую собеседуется ?
    salary_from         integer, -- нижняя граница зарплаты
    email               varchar(255), -- электронная почта
    telegram            varchar(255), -- телеграмм
    skype               varchar(255), -- skype, скайп
    cell_phone          varchar(255), -- телефон
    office_phone        varchar(255), -- рабочий телефон, редко заполняют, указывают второй телефон
    region_name         varchar(255), -- город
    country_name        varchar(255), -- страна
    created_at          timestamp(0), -- когда была создана заявка (формат временной)
    updated_at          timestamp(0), -- когда обновлялась заявка (формат временной)
    created_by          integer
        references "user", -- айдишник того, кто создал заявку, ссылается на табличку user
    updated_by          integer
        references "user", -- айдишник того, кто обновил заявку, ссылается на табличку user
    candidate_source_id integer not null
        references candidate_source
            on update cascade on delete cascade, -- 
    uuid                varchar(36), -- айдишник файла
    owned_at            timestamp(0), -- когда был закреплён?
    owned_by            integer -- айдишник того, чья это заявка, ссылается на табличку user
        references "user",
    last_moved_at       timestamp(0), -- в последний раз менялся ?
    last_moved_by       integer
        references "user", -- айдишник того, кто в последний раз менял, ссылается на табличку user
    chat_id             integer
        references chat -- айдишник комментария ?
);

alter table candidate
    owner to hrtools;

create index candidate_candidate_source_id_idx
    on candidate (candidate_source_id);

create index candidate_chat_id_idx
    on candidate (chat_id);

create index candidate_created_at_idx
    on candidate (created_at);

create index candidate_created_by_idx
    on candidate (created_by);

create index candidate_last_moved_by_idx
    on candidate (last_moved_by);

create index candidate_owned_at_idx
    on candidate (owned_at);

create index candidate_owned_by_idx
    on candidate (owned_by);

create index candidate_updated_by_idx
    on candidate (updated_by);

create unique index candidate_uuid_idx
    on candidate (uuid);

grant select on candidate to grafana;

grant select on candidate to clickhouse;


```



## 9. Таблица candidate_advanced_search
``` sql
-- auto-generated definition
create table candidate_advanced_search
(
    candidate_id        integer, -- айдишник кандидата
    experience_duration double precision, -- какой опыт (в чём считается опыт?)
    age                 double precision, -- возраст (в чём ?)
    last_moved_at       timestamp(0), -- в последний раз менял статус ?
    position_history    tsvector, -- кем успел поработать
    key_skills          tsvector, -- что умеет, поле может быть пустым
    description_summary tsvector -- ?
);

alter table candidate_advanced_search
    owner to hrtools;

create index candidate_advanced_search_age_idx
    on candidate_advanced_search (age);

create unique index candidate_advanced_search_candidate_id_idx
    on candidate_advanced_search (candidate_id);

create index candidate_advanced_search_description_summary_idx
    on candidate_advanced_search using gin (description_summary);

create index candidate_advanced_search_experience_duration_idx
    on candidate_advanced_search (experience_duration);

create index candidate_advanced_search_key_skills_idx
    on candidate_advanced_search using gin (key_skills);

create index candidate_advanced_search_last_moved_at_idx
    on candidate_advanced_search (last_moved_at);

create index candidate_advanced_search_position_history_idx
    on candidate_advanced_search using gin (position_history);

grant select on candidate_advanced_search to grafana;

```

## 10. Таблица candidate_comment

``` sql
-- auto-generated definition
create table candidate_comment
(
    candidate_comment_id serial
        primary key, -- айдишник комментария кандидата, первичный ключ, 
    comment              varchar(255), -- комментарий
    created_at           timestamp(0), -- когда был создан
    created_by           integer not null
        references "user", -- кем был создан, ссылка на айдишник в таблице user 
    candidate_id         integer not null
        references candidate 
            on update cascade on delete cascade -- айдишник кандидата, ссылается на айдишник в таблице candidate
);

alter table candidate_comment
    owner to hrtools;

create index candidate_comment_candidate_id_idx
    on candidate_comment (candidate_id);

create index candidate_comment_created_by_idx
    on candidate_comment (created_by);

grant select on candidate_comment to grafana;

grant select on candidate_comment to clickhouse;


```

## 11. Таблица candidate_contact_watch_list

``` sql

-- auto-generated definition
create table candidate_contact_watch_list
(
    candidate_contact_watch_list_id integer default nextval('candidate_contact_watch_list_candidate_contact_watch_list_i_seq'::regclass) not null
        primary key, -- первичный ключ, айдишник ?
    created_at                      timestamp(0), -- когда был создан
    user_id                         integer
        references "user", -- кем был создан (?), ссылается на айдишник в таблице user
    candidate_id                    integer                                                                                              not null
        references candidate
            on update cascade on delete cascade -- айдишник кандидата, ссылается на айдишник в таблице candidate
);

alter table candidate_contact_watch_list
    owner to hrtools;

create index candidate_contact_watch_list_candidate_id_idx
    on candidate_contact_watch_list (candidate_id);

create index candidate_contact_watch_list_user_id_idx
    on candidate_contact_watch_list (user_id);

grant select on candidate_contact_watch_list to grafana;


```

## 12. Таблица 

``` sql
-- auto-generated definition
create table candidate_file
(
    candidate_file_id      serial
        primary key, -- айдишник файла кандидата
    path                   varchar(255) not null, -- путь до файла
    created_at             timestamp(0), -- когда был создан
    created_by             integer
        references "user", -- кем был создан, ссылается на айдишник в табличке user
    filename               varchar(255), -- имя файла
    uuid                   varchar(36), -- айдишник файла ?
    candidate_id           integer
                                        references candidate
                                            on update cascade on delete set null, -- айдишник кандидата, ссылается на айдишник кандидата в таблице candidates
    updated_at             timestamp(0), -- когда было обновлено
    candidate_file_type_id integer      not null
        references candidate_file_type
            on update cascade on delete set null, -- тип файла, ссылается на табличку candidate_file_type
    updated_by             integer
        references "user" -- кем обновлялось, ссылка на айдишник пользователя в таблице users
);

alter table candidate_file
    owner to hrtools;

create index candidate_file_candidate_file_type_id_idx
    on candidate_file (candidate_file_type_id);

create index candidate_file_candidate_id_idx
    on candidate_file (candidate_id);

create index candidate_file_created_by_idx
    on candidate_file (created_by);

create index candidate_file_updated_by_idx
    on candidate_file (updated_by);

grant select on candidate_file to grafana;

grant select on candidate_file to clickhouse;


```

## 13. Таблица candidate_file_type

``` sql
-- auto-generated definition
create table candidate_file_type
(
    candidate_file_type_id serial
        primary key, -- айдишник файла кандидата
    name                   varchar(255), -- название файла, например "файл с резюме"
    shortcut               varchar(255) -- сокращение для названия файла
);

alter table candidate_file_type
    owner to hrtools;

grant select on candidate_file_type to grafana;

grant select on candidate_file_type to clickhouse;


```

## 14. Таблица candidate_responsible

``` sql
-- auto-generated definition
create table candidate_responsible
(
    candidate_responsible_id serial
        primary key, -- первичный ключ, айдишник (чей?)
    candidate_id             integer not null
        references candidate 
            on update cascade on delete cascade, -- айдишник кандидата, ссылается на табличку candidates
    user_id                  integer not null
        references "user"
            on update cascade on delete cascade -- айдишник пользователя, ссылается на табличку users
);

alter table candidate_responsible
    owner to hrtools;

create index candidate_responsible_candidate_id_idx
    on candidate_responsible (candidate_id);

create unique index candidate_responsible_candidate_id_user_id_idx
    on candidate_responsible (candidate_id, user_id);

create index candidate_responsible_user_id_idx
    on candidate_responsible (user_id);

grant select on candidate_responsible to grafana;

grant select on candidate_responsible to clickhouse;

```

## 15. Таблица candidate_source

``` sql
-- auto-generated definition
create table candidate_source
(
    candidate_source_id serial
        primary key, -- первичный ключ, айдишник источника кандидата
    name                varchar(255), -- название, откуда пришёел кандидат (например, Habr)
    description         varchar(255) -- описание каждого названия
);

alter table candidate_source
    owner to hrtools;

grant select on candidate_source to grafana;

grant select on candidate_source to clickhouse;



```

## 16. Таблица chat
``` sql
-- auto-generated definition
-- Видимо, это айдишники возможных комментов и их типов
create table chat
(
    chat_id      serial
        primary key, -- первичный ключ, айдишник комментария
    created_at   timestamp(0), -- когда был создан
    chat_type_id integer not null
        references chat_type -- тип комментария, ссылается на айдишник в табличке chat_type
);

alter table chat
    owner to hrtools;

create index chat_chat_type_id_idx
    on chat (chat_type_id);

grant select on chat to grafana;



```

## 17. Таблица chat_comment


``` sql
-- auto-generated definition
create table chat_comment
(
    chat_comment_id serial
        primary key, -- первичный ключ, айдишник комментария
    comment         varchar(255) not null, -- сам комментарий
    created_at      timestamp(0), -- когда был создан
    created_by      integer 
        references "user", -- кем был создан, ссылка на таблицу user
    chat_id         integer      not null
        references chat
);

alter table chat_comment
    owner to hrtools;

create index chat_comment_chat_id_idx
    on chat_comment (chat_id);

create index chat_comment_created_by_idx
    on chat_comment (created_by);

grant select on chat_comment to grafana;

```


## 17. Таблица chat_type


``` sql
-- auto-generated definition
create table chat_type
(
    chat_type_id serial
        primary key, -- айдишник коммента
    name         varchar(255) not null -- чей комментарий
);

alter table chat_type
    owner to hrtools;

grant select on chat_type to grafana;


```

## 18. Таблица client
``` sql
-- auto-generated definition
create table client
(
    client_id                    serial
        primary key, -- первичный ключ, айдишник клиента
    name                         varchar(255) not null, -- имя клиента
    shortcut                     varchar(255), -- сокращение имени
    address                      varchar(255), -- адрес (пока все пустые)
    comment                      varchar(255), -- комментарий
    created_at                   timestamp(0), -- когда был создан
    updated_at                   timestamp(0), -- когда был обновлён
    created_by                   integer
        references "user", -- кем был создан, ссылается на айдишник пользователя в табличке users
    client_status_id             integer      not null
        references client_status, -- айдишник статуса клиента, ссылается на табличку client_status
    logo_path                    varchar(255), -- путь до ?
    uuid                         varchar(36), -- айдишник файла ? 
    is_starred                   boolean default false, -- начат или нет ?
    default_pipeline_template_id integer
        references pipeline_template, -- айдишник этапа, ссылается на pipeline_template
    chat_id                      integer
        references chat -- айдишник комментарий
);

alter table client
    owner to hrtools;

create index client_chat_id_idx
    on client (chat_id);

create index client_client_status_id_idx
    on client (client_status_id);

create index client_created_by_idx
    on client (created_by);

create index client_default_pipeline_template_id_idx
    on client (default_pipeline_template_id);

grant select on client to grafana;

grant select on client to clickhouse;
```

## 19. Таблица client_project

``` sql
-- auto-generated definition
create table  client_project
(
    client_project_id    serial
        primary key, -- первичный ключ, айдишник проекта клиента
    name                 varchar(255), -- название проекта
    chief                varchar(255), -- имя начальника
    gross_amount         bigint, -- зарплата гросс
    start_date           timestamp(0), -- время начала
    end_date             timestamp(0), -- время конца
    created_at           timestamp(0), -- когда был создан 
    updated_at           timestamp(0), -- когда был обновлён
    created_by           integer
                                 references "user"
                                     on update cascade on delete set null, -- кем был создан, ссылается на айдишник пользователя в табличке users
    updated_by           integer
                                 references "user"
                                     on update cascade on delete set null, -- кем был обновлён, ссылается на айдишник пользователя в табличке users
    client_id            integer not null
        references client
            on update cascade on delete cascade, -- айдишник клиента, ссылается на айдишник в табличке clients
    payment_rate_type_id integer
                                 references payment_rate_type
                                     on update cascade on delete set null, -- айдишник типа оплаты ?, ссылается на айдишник в таблице payment_rate_type
    currency_id          integer
                                 references currency
                                     on update cascade on delete set null, -- айдишник валюты, ссылается на на айдишник в таблице currency
    user_id              integer
                                 references "user"
                                     on update cascade on delete set null -- айдишник пользователя, ссылается на айдишник в таблице users
);

alter table client_project
    owner to hrtools;

create index client_project_client_id_idx
    on client_project (client_id);

create index client_project_created_by_idx
    on client_project (created_by);

create index client_project_currency_id_idx
    on client_project (currency_id);

create index client_project_payment_rate_type_id_idx
    on client_project (payment_rate_type_id);

create index client_project_updated_by_idx
    on client_project (updated_by);

create index client_project_user_id_idx
    on client_project (user_id);

grant select on client_project to grafana;


```

## 20. Таблица client_status
``` sql
-- auto-generated definition
create table client_status
(
    client_status_id serial
        primary key, -- первичный ключ, айдишник статуса
    name             varchar(255) not null -- сам статус (клиент, архив и т.п.)
);

alter table client_status
    owner to hrtools;

grant select on client_status to grafana;

grant select on client_status to clickhouse;
```

## 21. Таблица client_stream

``` sql
-- auto-generated definition
create table client_stream
(
    client_stream_id serial
        primary key, -- первичный ключ, айдишник чего (?)
    name             varchar(255) not null, -- типы (например, 4К - КПиП ЦФТ, 4К - Кредитные продукты и процессы, Брокерское обслуживание юр.лиц)
    shortcut         varchar(255) not null, -- сокращение
    description      text, -- описание (не у всех есть)
    client_id        integer      not null 
        references client -- айдишник клиента, ссылается на айдишники в таблице client
);

alter table client_stream
    owner to hrtools;

create index client_stream_client_id_idx
    on client_stream (client_id);

grant select on client_stream to grafana;

grant select on client_stream to clickhouse;


```

## 22. Таблица contact_person
``` sql
-- auto-generated definition
create table contact_person
(
    contact_person_id        serial
        primary key, -- первичный ключ, айдишник контактного лица
    last_name                varchar(255), -- фамилия
    first_name               varchar(255)      not null, -- имя
    middle_name              varchar(255), -- отчество
    company                  varchar(255), -- название компании
    email                    varchar(255), -- электронная почта
    telegram                 varchar(255), -- телеграмм
    skype                    varchar(255), -- skype, скайп
    cell_phone               varchar(255), -- номер телефона
    office_phone             varchar(255), -- рабочий телефон (часто городской)
    office_phone_extension   varchar(255), -- рабочий телефон -- уточнение ?
    birth_date               date, -- дата рождения (пока пустая)
    comment                  varchar(255), -- комментарий
    created_at               timestamp(0), -- когда был создан
    updated_at               timestamp(0), -- когда был обновлён 
    created_by               integer
        references "user", -- кем был создан, ссылается на айдишник в табличке user
    client_id                integer
        references client, -- айдишник клиента, ссылается на айдишник в табличке client
    position                 varchar(255), -- должность
    contact_person_status_id integer default 1 not null
        references contact_person_status, -- айдишник статуса пользователя, ссылается на табличку contact_person_status
    chat_id                  integer
        references chat -- айдишник комментария, ссылается на табличку chat
);

alter table contact_person
    owner to hrtools;

create index contact_person_chat_id_idx
    on contact_person (chat_id);

create index contact_person_client_id_idx
    on contact_person (client_id);

create index contact_person_contact_person_status_id_idx
    on contact_person (contact_person_status_id);

create index contact_person_created_by_idx
    on contact_person (created_by);

grant select on contact_person to grafana;

grant select on contact_person to clickhouse;


```

## 23. Таблица contact_person_status

``` sql
-- auto-generated definition
create table contact_person_status
(
    contact_person_status_id serial
        primary key, -- первичный ключ, айдишник статуса
    name                     varchar(255) not null -- статус (работает, не работает, перешел к другому клиенту
);

alter table contact_person_status
    owner to hrtools;

grant select on contact_person_status to grafana;


```

## 24. Таблица contact_type
``` sql
-- auto-generated definition
create table contact_type
(
    contact_type_id serial
        primary key,-- первичный ключ, айдишник типа контакта
    name            varchar(255) not null -- тип контакта (skype, telegram)
);

alter table contact_type
    owner to hrtools;

grant select on contact_type to grafana;

grant select on contact_type to clickhouse;


```

## 25. Таблица country

``` sql
-- auto-generated definition
create table country
(
    country_id serial
        primary key, -- первичный ключ, айдишник страны
    name       varchar(255) not null -- название страны
);

alter table country
    owner to hrtools;

grant select on country to grafana;


```

## 26. Таблица cover_letter

``` sql
-- auto-generated definition
create table cover_letter
(
    cover_letter_id    serial
        primary key, -- первичный ключ, айдишник сопроводительного письма
    name               varchar(255), -- название сопровдительного письма
    experience         text, -- опыт
    skills             text, -- навыки
    motivation         text, -- мотивация
    additional         text, -- дополнительно
    city               text, -- город
    salary             text, -- зарплата (иногда с комментарием)
    interview_link     text, -- ссылка на интервью (все пустые)
    candidate_link     text, -- ссылка на кандидата ? (все пустые)
    meeting_option     text, -- тип встречи ?
    created_at         timestamp(0), -- когда были созданы
    updated_at         timestamp(0), -- когда были обновлены
    created_by         integer
        references "user", -- кем был создан, ссылается на айдишник в табличке users
    candidate_id       integer not null
        references candidate
            on update cascade on delete cascade, -- кто кандидат, ссылается на айдишник в табличке candidates
    updated_by         integer
        references "user", -- кем был обновлён, ссылается на айдишник в табличке users
    date_of_employment timestamp(0) -- дата трудоустройства (всё пусто)
);

alter table cover_letter
    owner to hrtools;

create index cover_letter_candidate_id_idx
    on cover_letter (candidate_id);

create index cover_letter_created_by_idx
    on cover_letter (created_by);

create index cover_letter_updated_by_idx
    on cover_letter (updated_by);

grant select on cover_letter to grafana;

grant select on cover_letter to clickhouse;


```
