# zahra_backend

npx sequelize-cli model:generate --name Request --attributes pick_up_date:date,comment:string,state:string,made_by:string,id_ingrediente:integer
npx sequelize-cli model:generate --name ShoppingCart --attributes email_cliente:string,id_ingredient:integer,estado:boolean
npx sequelize-cli model:generate --name Ingredient --attributes name:string, expiration_date:date, weight:integer,bought_date:date,price:float,description:string,owner:string
npx sequelize-cli model:generate --name User --attributes email:string,password:string,name:string,telephone:integer,member_since:date,address:string,description:string,is_admin:boolean
npx sequelize-cli model:generate --name Valoration --attributes comment:string,puntuation:integer,email_user:string,made_by:string



DOCUMENTATION API : 

AUTHENTIFICATION :

API Documentation
Endpoint /signup

HTTP Method: POST

Route: https://zahra-backend.onrender.com/signup

Arguments received: JSON object containing the following information:

email (string): User's email address.
password (string): User's password.
name (string): User's name.
telephone (string): User's telephone number.
member_since (string): User's membership start date.
address (string): User's address.
description (string): User's description.
is_admin (boolean): Flag indicating whether the user is an administrator or not.
Returns: JSON object with details of the created user:

email (string): Email address of the created user.
name (string): Name of the created user.
HTTP Status Codes:

201 Created: Successfully created a new user.
400 Bad Request: User with the specified email address already exists.
500 Internal Server Error: Server encountered an error while processing the request.
Endpoint /login

HTTP Method: POST

Route: https://zahra-backend.onrender.com/login

Arguments received: JSON object containing the following information:

email (string): User's email address for login.
password (string): User's password for login.
Returns: JSON object with details of successful login:

access_token (string): JWT token for user authentication.
token_type (string): Token type (Bearer).
expires_in (number): Number of seconds until the token expires.
user (object): Object containing details of the logged-in user:

email (string): Email address of the logged-in user.
name (string): Name of the logged-in user.
is_admin (boolean): Flag indicating whether the user is an administrator or not.
HTTP Status Codes:

200 OK: Successful login and JWT token returned.
401 Unauthorized: Incorrect password.
404 Not Found: User with the specified email address was not found.
500 Internal Server Error: Server encountered an error while processing the request.




API Documentation
Endpoint /login

HTTP Method: POST
Route: https://zahra-backend.onrender.com/login
Arguments Received: JSON object containing the following information:
email (string): The user's email address for login.
password (string): The user's password for login.
Returns: JSON object with details of successful login:
access_token (string): JWT token for user authentication.
token_type (string): Type of token (Bearer).
expires_in (number): Number of seconds until the token expires.
user (object): Object containing details of the logged-in user:
email (string): The user's email address.
name (string): The user's name.
is_admin (boolean): Indicator if the user is an administrator or not.
HTTP Status Codes:
200 OK: Successful login and JWT token returned.
401 Unauthorized: Incorrect password.
404 Not Found: The user with the specified email address was not found.
500 Internal Server Error: Server error occurred while processing the request.



_______



INGREDIENTES


API Documentation
Endpoint /ingredientes

Create Ingredient

HTTP Method: POST
Route: https://zahra-backend.onrender.com/ingredientes
Arguments Received: JSON object containing the following information:
name (string): The name of the ingredient.
description (string): The description of the ingredient.
owner (string): The owner or user associated with the ingredient.
Returns: JSON object with details of the created ingredient:
id (number): The unique identifier of the ingredient.
name (string): The name of the ingredient.
description (string): The description of the ingredient.
owner (string): The owner or user associated with the ingredient.
HTTP Status Codes:
201 Created: Successfully created a new ingredient.
500 Internal Server Error: Server error occurred while processing the request.
Delete Ingredient by ID

HTTP Method: DELETE
Route: https://zahra-backend.onrender.com/ingredientes/:id
Arguments Received: Path parameter id (number) indicating the ID of the ingredient to delete.
Returns: No content.
HTTP Status Codes:
204 No Content: Successfully deleted the ingredient.
404 Not Found: The ingredient with the specified ID was not found.
500 Internal Server Error: Server error occurred while processing the request.
Get All Ingredients

HTTP Method: GET
Route: https://zahra-backend.onrender.com/ingredientes
Returns: JSON array containing all ingredients.
HTTP Status Codes:
200 OK: Successfully retrieved all ingredients.
500 Internal Server Error: Server error occurred while processing the request.
Update Ingredient by ID

HTTP Method: PATCH
Route: https://zahra-backend.onrender.com/ingredientes/:id
Arguments Received: Path parameter id (number) indicating the ID of the ingredient to update. JSON object containing fields to update.
Returns: JSON object with details of the updated ingredient.
HTTP Status Codes:
200 OK: Successfully updated the ingredient.
404 Not Found: The ingredient with the specified ID was not found.
500 Internal Server Error: Server error occurred while processing the request.
Get Ingredients by Owner's Email

HTTP Method: GET
Route: https://zahra-backend.onrender.com/ingredientes/user/:email
Arguments Received: Path parameter email (string) indicating the email address of the ingredient owner.
Returns: JSON array containing all ingredients associated with the specified owner's email.
HTTP Status Codes:
200 OK: Successfully retrieved ingredients by owner's email.
500 Internal Server Error: Server error occurred while processing the request.
Get Ingredient by ID

HTTP Method: GET
Route: https://zahra-backend.onrender.com/ingredientes/:id
Arguments Received: Path parameter id (number) indicating the ID of the ingredient to retrieve.
Returns: JSON object with details of the specified ingredient.
HTTP Status Codes:
200 OK: Successfully retrieved the ingredient.
404 Not Found: The ingredient with the specified ID was not found.
500 Internal Server Error: Server error occurred while processing the request.



____________

REQUESTS



API Documentation
Endpoint /requests

Create Request

HTTP Method: POST
Route: https://zahra-backend.onrender.com/requests
Arguments Received: JSON object containing the following information:
pick_up_date (string): The date for picking up the request.
comment (string): Optional comment associated with the request.
state (string): State of the request.
made_by (string): Email of the user making the request.
id_ingrediente (number): ID of the ingredient related to the request.
Returns: JSON object with details of the created request:
id (number): Unique identifier of the request.
pick_up_date (string): The date for picking up the request.
comment (string): Optional comment associated with the request.
state (string): State of the request.
made_by (string): Email of the user making the request.
id_ingrediente (number): ID of the ingredient related to the request.
HTTP Status Codes:
201 Created: Successfully created a new request.
500 Internal Server Error: Server error occurred while processing the request.
Get Requests by Ingredient ID

HTTP Method: GET
Route: https://zahra-backend.onrender.com/requests/ingredient/:id_ingrediente
Arguments Received: Path parameter id_ingrediente (number) indicating the ID of the ingredient to retrieve requests for.
Returns: JSON array containing all requests associated with the specified ingredient ID.
HTTP Status Codes:
200 OK: Successfully retrieved requests by ingredient ID.
500 Internal Server Error: Server error occurred while processing the request.
Update Request by Ingredient ID and Request ID

HTTP Method: PATCH
Route: https://zahra-backend.onrender.com/requests/ingredient/:id_ingrediente
Arguments Received: Path parameter id_ingrediente (number) indicating the ID of the ingredient related to the request. JSON object containing fields to update:
id (number): ID of the request to update.
pick_up_date (string): The date for picking up the request.
comment (string): Optional comment associated with the request.
state (string): State of the request.
made_by (string): Email of the user making the request.
Returns: JSON object with details of the updated request.
HTTP Status Codes:
200 OK: Successfully updated the request.
404 Not Found: The request with the specified IDs was not found.
500 Internal Server Error: Server error occurred while processing the request.
Get Requests by User Email

HTTP Method: GET
Route: https://zahra-backend.onrender.com/requests/user/:made_by
Arguments Received: Path parameter made_by (string) indicating the email address of the user to retrieve requests for.
Returns: JSON array containing all requests associated with the specified user email.
HTTP Status Codes:
200 OK: Successfully retrieved requests by user email.
500 Internal Server Error: Server error occurred while processing the request.
Delete Request by User Email and Ingredient ID

HTTP Method: DELETE
Route: https://zahra-backend.onrender.com/requests/user/:made_by/ingredient/:id_ingrediente
Arguments Received: Path parameters made_by (string) indicating the email address of the user and id_ingrediente (number) indicating the ID of the ingredient related to the request.
Returns: No content.
HTTP Status Codes:
204 No Content: Successfully deleted the request.
404 Not Found: The request with the specified user email and ingredient ID was not found.
500 Internal Server Error: Server error occurred while processing the request.


________

USERS



API Documentation
Endpoint /users

Create User

HTTP Method: POST
Route: https://zahra-backend.onrender.com/users
Arguments Received: JSON object containing the following information:
email (string): Email address of the user.
password (string): Password of the user.
name (string): Name of the user.
telephone (string): Telephone number of the user.
member_since (string): Date of membership of the user.
address (string): Address of the user.
description (string): Description of the user.
is_admin (boolean): Flag indicating if the user is an administrator.
Returns: JSON object with details of the created user:
id (number): Unique identifier of the user.
email (string): Email address of the user.
password (string): Password hash of the user.
name (string): Name of the user.
telephone (string): Telephone number of the user.
member_since (string): Date of membership of the user.
address (string): Address of the user.
description (string): Description of the user.
is_admin (boolean): Flag indicating if the user is an administrator.
HTTP Status Codes:
201 Created: Successfully created a new user.
400 Bad Request: Client error occurred due to invalid input data.
500 Internal Server Error: Server error occurred while processing the request.
Get All Users

HTTP Method: GET
Route: https://zahra-backend.onrender.com/users
Returns: JSON array containing all users stored in the database.
HTTP Status Codes:
200 OK: Successfully retrieved all users.
500 Internal Server Error: Server error occurred while processing the request.
Get User by Email

HTTP Method: GET
Route: https://zahra-backend.onrender.com/users/:email
Arguments Received: Path parameter email (string) indicating the email address of the user to retrieve.
Returns: JSON object with details of the user matching the provided email:
id (number): Unique identifier of the user.
email (string): Email address of the user.
password (string): Password hash of the user.
name (string): Name of the user.
telephone (string): Telephone number of the user.
member_since (string): Date of membership of the user.
address (string): Address of the user.
description (string): Description of the user.
is_admin (boolean): Flag indicating if the user is an administrator.
HTTP Status Codes:
200 OK: Successfully retrieved the user by email.
404 Not Found: User with the specified email address was not found.
500 Internal Server Error: Server error occurred while processing the request.
Delete User by Email

HTTP Method: DELETE
Route: https://zahra-backend.onrender.com/users/deleteByEmail/:email
Arguments Received: Path parameter email (string) indicating the email address of the user to delete.
Returns: No content.
HTTP Status Codes:
204 No Content: Successfully deleted the user.
404 Not Found: User with the specified email address was not found.
500 Internal Server Error: Server error occurred while processing the request.
Update User by ID

HTTP Method: PATCH
Route: https://zahra-backend.onrender.com/users/:id
Arguments Received: Path parameter id (number) indicating the ID of the user to update. JSON object containing fields to update:
email (string): Email address of the user.
password (string): Password hash of the user.
name (string): Name of the user.
telephone (string): Telephone number of the user.
member_since (string): Date of membership of the user.
address (string): Address of the user.
description (string): Description of the user.
is_admin (boolean): Flag indicating if the user is an administrator.
Returns: JSON object with details of the updated user.
HTTP Status Codes:
200 OK: Successfully updated the user.
404 Not Found: User with the specified ID was not found.
500 Internal Server Error: Server error occurred while processing the request.



______________


VALORATIONS



API Documentation
Endpoint /valorations

Create Valoration

HTTP Method: POST
Route: https://zahra-backend.onrender.com/valorations
Arguments Received: JSON object containing the following information:
comment (string): Comment associated with the valoration.
puntuation (number): Rating given by the user (1 to 5).
email_user (string): Email address of the user receiving the valoration.
made_by (string): Email address of the user who made the valoration.
Returns: JSON object with details of the created valoration:
id (number): Unique identifier of the valoration.
comment (string): Comment associated with the valoration.
puntuation (number): Rating given by the user (1 to 5).
email_user (string): Email address of the user receiving the valoration.
made_by (string): Email address of the user who made the valoration.
HTTP Status Codes:
201 Created: Successfully created a new valoration.
500 Internal Server Error: Server error occurred while processing the request.
Get Valoration by Email User

HTTP Method: GET
Route: https://zahra-backend.onrender.com/valorations/:email_user
Arguments Received: Path parameter email_user (string) indicating the email address of the user to retrieve valuations for.
Returns: JSON array containing all valuations made for the specified user:
id (number): Unique identifier of the valoration.
comment (string): Comment associated with the valoration.
puntuation (number): Rating given by the user (1 to 5).
email_user (string): Email address of the user receiving the valoration.
made_by (string): Email address of the user who made the valoration.
HTTP Status Codes:
200 OK: Successfully retrieved valuations for the specified user.
404 Not Found: No valuations found for the specified user.
500 Internal Server Error: Server error occurred while processing the request.
Update Valoration by ID

HTTP Method: PATCH
Route: https://zahra-backend.onrender.com/valorations/:id
Arguments Received: Path parameter id (number) indicating the ID of the valoration to update. JSON object containing fields to update:
comment (string): Comment associated with the valoration.
puntuation (number): Rating given by the user (1 to 5).
email_user (string): Email address of the user receiving the valoration.
made_by (string): Email address of the user who made the valoration.
Returns: JSON object with details of the updated valoration.
HTTP Status Codes:
200 OK: Successfully updated the valoration.
404 Not Found: Valoration with the specified ID was not found.
500 Internal Server Error: Server error occurred while processing the request.
Delete Valoration by User

HTTP Method: DELETE
Route: https://zahra-backend.onrender.com/valorations/user/:made_by
Arguments Received: Path parameter made_by (string) indicating the email address of the user who made the valoration to delete.
Returns: No content.
HTTP Status Codes:
204 No Content: Successfully deleted the valoration.
404 Not Found: Valoration made by the specified user was not found.
500 Internal Server Error: Server error occurred while processing the request.
Check Existing Valoration

HTTP Method: GET
Route: https://zahra-backend.onrender.com/valorations
Query Parameters: email_user (string) and made_by (string) to check if a valoration exists between these two users.
Returns: JSON object with details of the existing valoration if found:
id (number): Unique identifier of the valoration.
comment (string): Comment associated with the valoration.
puntuation (number): Rating given by the user (1 to 5).
email_user (string): Email address of the user receiving the valoration.
made_by (string): Email address of the user who made the valoration.
HTTP Status Codes:
200 OK: Successfully found an existing valoration.
404 Not Found: No valoration found between the specified users.
500 Internal Server Error: Server error occurred while processing the request.











