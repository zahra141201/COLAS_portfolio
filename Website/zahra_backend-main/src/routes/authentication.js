const Router = require('koa-router');
const jwt = require('jsonwebtoken');
const bcrypt = require('bcrypt');
const dotenv = require('dotenv');
const { User } = require('../models'); // Assurez-vous d'importer correctement votre modèle User

dotenv.config();

const router = new Router();

router.post("/signup", async (ctx) => {
    const authInfo = ctx.request.body;
    console.log('Received authInfo:', authInfo);
    try {
        ctx.body= 'coucou'
        ctx.body = `The user has email '${authInfo.email}'`
        let user = await User.findOne({ where: { email: authInfo.email } });
        if (user) {
            ctx.body = `The user with email '${authInfo.email}' already exists`;
            ctx.status = 400;
            return;
        }

        const saltRounds = 10;
        const hashPassword = await bcrypt.hash(authInfo.password, saltRounds);

        user = await User.create({
            email: authInfo.email,
            password: hashPassword,
            name: authInfo.name,
            telephone: authInfo.telephone,
            member_since: authInfo.member_since,
            address: authInfo.address,
            description: authInfo.description,
            is_admin: authInfo.is_admin
        });

        ctx.body = {
            email: user.email,
            name: user.name
        };
        ctx.status = 201;
    } catch (error) {
        console.error('Error creating user:', error);
        ctx.body = error.message || 'An error occurred while processing your request';
        ctx.status = 500;
    }
});

router.post("/login", async (ctx) => {
    const authInfo = ctx.request.body;
    try {
        let user = await User.findOne({ where: { email: authInfo.email } });
        if (!user) {
            ctx.body = `The user with email '${authInfo.email}' was not found`;
            ctx.status = 404; // Utilisez le code approprié pour l'utilisateur non trouvé
            return;
        }

        const validPassword = await bcrypt.compare(authInfo.password, user.password);
        if (validPassword) {
            const expirationSeconds = 1 * 60 * 60 * 24;
            const JWT_PRIVATE_KEY = process.env.JWT_SECRET;

            const token = jwt.sign(
                { scope: ['user'] ,
                is_admin: user.is_admin,
        
                },
                JWT_PRIVATE_KEY,
                { subject: user.id.toString() },
                { expiresIn: expirationSeconds }
            );

            ctx.body = {
                access_token: token,
                token_type: 'Bearer',
                expires_in: expirationSeconds,
                user: {
                    email: user.email,
                    name: user.name,
                    is_admin: user.is_admin
                }
            };
            ctx.status = 200;
        } else {
            ctx.body = 'Incorrect password';
            ctx.status = 401; // Utilisez le code approprié pour un accès non autorisé
        }
    } catch (error) {
        console.error('Error logging in:', error);
        ctx.body = error.message || 'An error occurred while processing your request';
        ctx.status = 500;
    }
});

module.exports = router;
