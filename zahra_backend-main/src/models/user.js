'use strict';
const { Model } = require('sequelize');

module.exports = (sequelize, DataTypes) => {
  class User extends Model {
    static associate(models) {
      this.hasMany(models.Ingredient, {
        foreignKey: 'owner',
        onDelete: 'CASCADE',
        onUpdate: 'CASCADE', // Supprime les ingrédients associés lorsque l'utilisateur est supprimé
      });
      this.hasMany(models.Request, {
        foreignKey: 'made_by',
        onDelete: 'CASCADE',
        onUpdate: 'CASCADE', // Supprime les demandes associées lorsque l'utilisateur est supprimé
      });
      this.hasMany(models.Valoration, {
        foreignKey: 'email_user',
        onDelete: 'CASCADE',
        onUpdate: 'CASCADE', // Supprime les évaluations associées lorsque l'utilisateur est supprimé
      });
      this.hasMany(models.Valoration, {
        foreignKey: 'made_by',
        onDelete: 'CASCADE',
        onUpdate: 'CASCADE', // Supprime les évaluations faites par l'utilisateur lorsque l'utilisateur est supprimé
      });
      this.hasMany(models.ShoppingCart, {
        foreignKey: 'email_cliente',
        onDelete: 'CASCADE',
        onUpdate: 'CASCADE', // Supprime les paniers associés lorsque l'utilisateur est supprimé
      });
    }
  }

  User.init({
    email: DataTypes.STRING,
    password: DataTypes.STRING,
    name: DataTypes.STRING,
    telephone: DataTypes.INTEGER,
    member_since: DataTypes.DATE,
    address: DataTypes.STRING,
    description: DataTypes.STRING,
    is_admin: DataTypes.BOOLEAN
  }, {
    sequelize,
    modelName: 'User',
  });

  return User;
};
