'use strict';
const { Model } = require('sequelize');

module.exports = (sequelize, DataTypes) => {
  class Ingredient extends Model {
    static associate(models) {
      this.belongsTo(models.User, {
        foreignKey: 'owner',
        onDelete: 'CASCADE',
        onUpdate: 'CASCADE', // Supprime les ingrédients associés lorsque l'utilisateur est supprimé
      });
      this.hasMany(models.Request, {
        foreignKey: 'id_ingrediente',
        onDelete: 'CASCADE',
        onUpdate: 'CASCADE', // Supprime les demandes associées lorsque l'ingrédient est supprimé
      });
      this.hasMany(models.ShoppingCart, {
        foreignKey: 'id_ingredient',
        onDelete: 'CASCADE',
        onUpdate: 'CASCADE', // Supprime les entrées de panier associées lorsque l'ingrédient est supprimé
      });
    }
  }

  Ingredient.init({
    name: DataTypes.STRING,
    expiration_date: DataTypes.DATE,
    weight: DataTypes.INTEGER,
    bought_date: DataTypes.DATE,
    price: DataTypes.FLOAT,
    description: DataTypes.STRING,
    owner: DataTypes.STRING
  }, {
    sequelize,
    modelName: 'Ingredient',
  });

  return Ingredient;
};
